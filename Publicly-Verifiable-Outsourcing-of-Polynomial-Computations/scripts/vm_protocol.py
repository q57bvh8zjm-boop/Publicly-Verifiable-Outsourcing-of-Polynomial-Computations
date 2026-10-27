"""Framing, deployment metadata, and public experiment fixtures for TCP runs."""
import hashlib
import json
import os
from pathlib import Path
import platform
import socket
import struct
import subprocess

ROOT = Path(__file__).resolve().parents[1]
VERSION = 1
LIMIT = 512 * 1024 * 1024
SETUP, READY, QUERY, LOADED, GO, RESULT, BYE = range(1, 8)
BATCH_SETUP, BATCH_READY, BATCH_QUERY, BATCH_RESULT, ROW, DONE = range(20, 26)
HELLO, HELLO_REPLY, ERROR = 100, 101, 102
U32, U64, HEADER = struct.Struct('!I'), struct.Struct('!Q'), struct.Struct('!IQ')
SCHEMES = ('Gamma_FE_Sh', 'Gamma_OV_Sh', 'Gamma_OV_HSS', 'Pi3')


def exact(stream, n):
    chunks = bytearray()
    while len(chunks) < n:
        part = stream.recv(n - len(chunks)) if isinstance(stream, socket.socket) else stream.read(n - len(chunks))
        if not part:
            raise EOFError('connection or child process closed before a complete frame')
        chunks.extend(part)
    return bytes(chunks)


def recv_frame(stream, expected=None):
    kind, size = HEADER.unpack(exact(stream, HEADER.size))
    if size > LIMIT:
        raise ValueError('frame exceeds 512 MiB limit')
    payload = exact(stream, size)
    if kind == ERROR:
        raise RuntimeError(payload.decode('utf8', 'replace'))
    if expected is not None and kind != expected:
        raise ValueError(f'expected frame {expected}, received {kind}')
    return kind, payload


def send_frame(stream, kind, payload=b''):
    if len(payload) > LIMIT:
        raise ValueError('frame exceeds 512 MiB limit')
    header = HEADER.pack(kind, len(payload))
    if isinstance(stream, socket.socket):
        stream.sendall(header)
        stream.sendall(payload)
    else:
        stream.write(header)
        stream.write(payload)
        stream.flush()


def unpack_batch(data):
    if len(data) < 4:
        raise ValueError('short batch')
    count = U32.unpack_from(data)[0]
    if not 1 <= count <= 128:
        raise ValueError('invalid server count')
    pos, items = 4, []
    for expected in range(count):
        if len(data) - pos < 12:
            raise ValueError('short batch entry')
        index, size = HEADER.unpack_from(data, pos)
        pos += 12
        if index != expected or size > len(data) - pos:
            raise ValueError('invalid batch order or length')
        items.append(data[pos:pos + size])
        pos += size
    if pos != len(data):
        raise ValueError('trailing batch bytes')
    return items


def pack_batch(items):
    return U32.pack(len(items)) + b''.join(HEADER.pack(i, len(p)) + p for i, p in enumerate(items))


def source_hash():
    h = hashlib.sha256()
    paths = [ROOT / 'Makefile']
    for folder in ('lib', 'protocols', 'benchmarks', 'scripts', 'tests'):
        paths += [p for p in (ROOT / folder).rglob('*') if p.suffix in ('.c', '.h', '.py', '.sh')]
    for p in sorted(paths):
        h.update(str(p.relative_to(ROOT)).encode() + b'\0' + p.read_bytes() + b'\0')
    return h.hexdigest()


def command(*args):
    try:
        p = subprocess.run(args, cwd=ROOT, text=True, capture_output=True, timeout=10)
        return p.stdout.strip() if p.returncode == 0 else {'unavailable': p.stderr.strip() or p.returncode}
    except (OSError, subprocess.TimeoutExpired) as e:
        return {'unavailable': str(e)}


def read_optional(path):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def metadata():
    frequencies = {}
    for p in sorted(Path('/sys/devices/system/cpu').glob('cpu[0-9]*/cpufreq/*')):
        if p.name in ('scaling_governor', 'scaling_cur_freq', 'scaling_min_freq', 'scaling_max_freq'):
            frequencies[str(p)] = read_optional(p)
    return {
        'wire_version': VERSION, 'source_sha256': source_hash(),
        'build': json.loads((ROOT / 'build/experiment-build.json').read_text()) if (ROOT / 'build/experiment-build.json').exists() else None,
        'hostname': platform.node(), 'platform': platform.platform(),
        'machine': platform.machine(), 'python': platform.python_version(),
        'logical_cpus': os.cpu_count(),
        'affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
        'cpuinfo': read_optional('/proc/cpuinfo'), 'meminfo': read_optional('/proc/meminfo'),
        'cpu_frequency': frequencies,
        'cgroup_cpu_max': read_optional('/sys/fs/cgroup/cpu.max'),
        'cgroup_memory_max': read_optional('/sys/fs/cgroup/memory.max'),
        'mac_cpu': command('sysctl', '-n', 'machdep.cpu.brand_string') if platform.system() == 'Darwin' else None,
        'mac_memory': command('sysctl', '-n', 'hw.memsize') if platform.system() == 'Darwin' else None,
        'git_revision': command('git', 'rev-parse', 'HEAD'),
        'compiler': command(os.environ.get('CC', 'cc'), '--version'),
        'build_commands': command('make', '-Bn', 'experiments'),
        'libraries': command('pkg-config', '--modversion', 'flint', 'gmp', 'libsodium'),
    }


def verify_build(binary):
    path = ROOT / 'build/experiment-build.json'
    if not path.exists():
        raise ValueError('missing build manifest; run make experiments')
    build = json.loads(path.read_text())
    actual = hashlib.sha256((ROOT / 'build' / binary).read_bytes()).hexdigest()
    if build['source_sha256'] != source_hash() or build['binaries'].get(binary) != actual:
        raise ValueError('source or binary changed since build; run make experiments')


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temp.replace(path)
