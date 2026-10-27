#!/usr/bin/env python3
"""Local TCP integration, malformed frames, identity rejection, and resume checks."""
import argparse
import csv
import json
import shutil
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from vm_protocol import HEADER, LIMIT, HELLO, HELLO_REPLY, VERSION, SETUP, recv_frame, send_frame


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--hss', action='store_true')
    ap.add_argument('--keep-results', type=Path)
    args = ap.parse_args()
    processes, logs = [], []
    with tempfile.TemporaryDirectory(prefix='pv-tcp-test-') as name:
        folder = Path(name)
        ports = []
        try:
            for i in range(3):
                with socket.socket() as sock:
                    sock.bind(('127.0.0.1', 0))
                    ports.append(sock.getsockname()[1])
                log = (folder / f'server{i}.log').open('w')
                logs.append(log)
                p = subprocess.Popen([sys.executable, str(ROOT / 'scripts/vm_server.py'),
                                      '--server-id', str(i), '--port', str(ports[i]), '--job-timeout', '180'],
                                     stdout=log, stderr=log)
                processes.append(p)
                for attempt in range(100):
                    try:
                        with socket.create_connection(('127.0.0.1', ports[i]), timeout=0.1):
                            break
                    except OSError:
                        if p.poll() is not None:
                            raise RuntimeError('test server exited')
                        time.sleep(0.05)
                else:
                    raise RuntimeError('server startup timeout')
            # The configured identity is enforced before starting cryptographic work.
            with socket.create_connection(('127.0.0.1', ports[0])) as conn:
                send_frame(conn, HELLO, json.dumps({'wire_version': VERSION, 'server_id': 1}).encode())
                try:
                    recv_frame(conn)
                    raise AssertionError('wrong server identity accepted')
                except RuntimeError as e:
                    assert 'identity' in str(e)
            # Oversized frames and malformed SETUP terminate only that C worker.
            for malformed in (HEADER.pack(SETUP, LIMIT + 1), HEADER.pack(SETUP, 1) + b'\0'):
                with socket.create_connection(('127.0.0.1', ports[0])) as conn:
                    conn.settimeout(15)
                    send_frame(conn, HELLO, json.dumps({'wire_version': VERSION, 'server_id': 0}).encode())
                    recv_frame(conn, HELLO_REPLY)
                    conn.sendall(malformed)
                    try:
                        recv_frame(conn)
                        raise AssertionError('malformed SETUP accepted')
                    except (EOFError, ConnectionResetError):
                        pass
                time.sleep(0.1)
            config = {'environment_label': 'local-tcp-functional-test', 'require_distinct_hosts': False,
                      'require_client_separate': False,
                      'servers': [{'server_id': i, 'host': '127.0.0.1', 'port': p} for i, p in enumerate(ports)]}
            (folder / 'config.json').write_text(json.dumps(config))
            schemes = ['Gamma_FE_Sh', 'Gamma_OV_Sh', 'Pi3']
            if args.hss:
                schemes.append('Gamma_OV_HSS')
            command = [sys.executable, str(ROOT / 'scripts/run_vm_experiments.py'), '--config', str(folder / 'config.json'),
                       '--profile', 'smoke', '--trials', '2', '--warmups', '0', '--native-repeats', '5',
                       '--timeout', '180', '--output', str(folder / 'results'), '--schemes', *schemes]
            subprocess.run(command, check=True)
            samples = []
            for path in (folder / 'results/jobs').glob('*/raw.csv'):
                with path.open() as f:
                    samples += list(csv.DictReader(f))
            assert len(samples) == 2 * len(schemes)
            for row in samples:
                times = json.loads(row['server_ms'])
                assert abs(float(row['server_max_ms']) - max(times)) < 1e-6
                assert abs(float(row['server_sum_ms']) - sum(times)) < 1e-6
                assert float(row['request_wire_bytes']) > 0 and float(row['response_wire_bytes']) > 0
                assert row['correct'] == '1'
                assert abs(float(row['client_ms']) - sum(float(row[k]) for k in ('probgen_ms', 'verify_ms', 'decode_ms'))) < 1e-6
            before = {str(p): p.read_bytes() for p in (folder / 'results/jobs').glob('*/raw.csv')}
            subprocess.run(command + ['--resume'], check=True)
            assert all(Path(p).read_bytes() == value for p, value in before.items())
            if args.keep_results:
                shutil.copytree(folder / 'results', args.keep_results)
            print(f'TCP integration passed: {len(samples)} correct samples; malformed/identity rejection; Ts=max; resume.')
        finally:
            for p in processes:
                p.terminate()
            for p in processes:
                p.wait(timeout=10)
            for log in logs:
                log.close()


if __name__ == '__main__':
    main()
