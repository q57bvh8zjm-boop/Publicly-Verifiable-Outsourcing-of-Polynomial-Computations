#!/usr/bin/env python3
"""Run real C protocols on TCP-connected VMs; preserve every measured sample."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
from pathlib import Path
import socket
import subprocess
import threading
import time
from vm_protocol import (ROOT, VERSION, SCHEMES, HELLO, HELLO_REPLY, SETUP, READY, QUERY,
                         LOADED, GO, RESULT, BYE, BATCH_SETUP, BATCH_READY, BATCH_QUERY,
                         BATCH_RESULT, ROW, DONE, U32, U64, HEADER, metadata, source_hash,
                         recv_frame, send_frame, unpack_batch, pack_batch, save_json, verify_build)


def jobs_for(profile, schemes):
    jobs = set()

    def add(names, ms, ds, ts):
        for s in names:
            for m in ms:
                for d in ds:
                    for t in ts:
                        jobs.add((s, m, d, t))

    if profile == 'smoke':
        add(SCHEMES, [2], [2], [1])
    if profile in ('figure2', 'paper'):
        add(SCHEMES[:3], range(200, 2001, 200), [2], [1])
    if profile in ('figure3', 'paper'):
        add(SCHEMES[:2], range(200, 2001, 200), [2], [2, 3])
    if profile in ('figure4', 'paper'):
        add([SCHEMES[0], 'Pi3'], range(200, 2001, 200), [2, 3], [1, 2])
    if profile in ('figure5', 'paper'):
        add(SCHEMES[:3], range(1, 16), [4, 8, 10], [1])
    return sorted(j for j in jobs if j[0] in schemes)


def identifier(job):
    s, m, d, t = job
    return f'{s}-m{m}-d{d}-t{t}'


def plan_job(job, limit, endpoints):
    s, m, d, t = job
    count = math.comb(m + d, d)
    servers = 2 if s == 'Gamma_OV_HSS' else d * t + 1
    reasons = []
    if count > limit:
        reasons.append(f'{count} dense coefficients exceeds limit {limit}')
    if servers > len(endpoints):
        reasons.append(f'needs {servers} servers; only {len(endpoints)} configured')
    # A field coefficient alone requires 32 bytes. FLINT containers, GMP limbs,
    # traversal temporaries, inputs and responses add to that lower bound.
    return {'id': identifier(job), 'scheme': s, 'm': m, 'd': d, 't': t,
            'servers': servers, 'coefficients': count,
            'coefficient_payload_min_bytes_per_node': count * 32,
            'planning_estimate_bytes_per_node': count * 128,
            'planning_estimate_is_not_a_memory_guarantee': True,
            'blocked': reasons}


def write_rows(path, rows):
    if not rows:
        return
    with Path(path).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, list) else v for k, v in row.items()})


def run_job(job, args, config, source, directory):
    s, m, d, t = job
    count = 2 if s == 'Gamma_OV_HSS' else d * t + 1
    endpoints = config['servers'][:count]
    sockets, rows, nodes = [], [], []
    process, timer = None, None
    expired = threading.Event()
    directory.mkdir(parents=True, exist_ok=True)
    command = [str(ROOT / 'build/vm_client'), s, str(m), str(d), str(t), str(args.trials),
               str(args.warmups), str(args.native_repeats), str(args.max_coefficients)]
    status = {'status': 'running', 'command': command, 'started': time.time(), 'source_sha256': source}
    save_json(directory / 'status.json', status)

    def expire():
        expired.set()
        if process is not None and process.poll() is None:
            process.kill()
        for conn in sockets:
            try:
                conn.shutdown(socket.SHUT_RDWR)
            except OSError:
                pass

    try:
        for i, endpoint in enumerate(endpoints):
            conn = socket.create_connection((endpoint['host'], endpoint['port']), timeout=30)
            sockets.append(conn)
            conn.settimeout(args.timeout)
            conn.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            send_frame(conn, HELLO, json.dumps({'wire_version': VERSION, 'server_id': i}).encode())
            _, response = recv_frame(conn, HELLO_REPLY)
            node = json.loads(response)
            if node.get('source_sha256') != source or node.get('server_id') != i:
                raise ValueError(f'server {i}: source or identity differs from client; deploy the same source archive')
            if node['max_coefficients'] < math.comb(m + d, d):
                raise ValueError(f'server {i}: coefficient limit is too small')
            nodes.append(node)
        save_json(directory / 'servers.json', nodes)
        if config.get('require_distinct_hosts', True):
            names = [n['hostname'] for n in nodes]
            if len(set(names)) != len(names):
                raise ValueError('server hostnames are not distinct; deploy one server on each VM')
            if config.get('require_client_separate', True) and socket.gethostname() in names:
                raise ValueError('client and server hostnames coincide; use a separate client VM')
        with (directory / 'client.stderr').open('w') as err, ThreadPoolExecutor(max_workers=count) as pool:
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err)
            timer = threading.Timer(args.timeout, expire)
            timer.daemon = True
            timer.start()

            def exchange(index, kind, payload, reply):
                try:
                    send_frame(sockets[index], kind, payload)
                    return recv_frame(sockets[index], reply)[1]
                except Exception:
                    # Wake every pending socket operation before joining the pool.
                    for conn in sockets:
                        try:
                            conn.shutdown(socket.SHUT_RDWR)
                        except OSError:
                            pass
                    raise

            while True:
                kind, payload = recv_frame(process.stdout)
                if kind == BATCH_SETUP:
                    parts = unpack_batch(payload)
                    if len(parts) != count:
                        raise ValueError('client/server count mismatch')
                    begin = time.perf_counter_ns()
                    futures = [pool.submit(exchange, i, SETUP, p, READY) for i, p in enumerate(parts)]
                    replies = [f.result() for f in futures]
                    elapsed = time.perf_counter_ns() - begin
                    send_frame(process.stdin, BATCH_READY, U64.pack(elapsed) + pack_batch(replies))
                elif kind == BATCH_QUERY:
                    parts = unpack_batch(payload)
                    if len(parts) != count:
                        raise ValueError('query count mismatch')
                    begin = time.perf_counter_ns()
                    futures = [pool.submit(exchange, i, QUERY, p, LOADED) for i, p in enumerate(parts)]
                    loaded = [f.result() for f in futures]
                    upload = time.perf_counter_ns() - begin
                    if any(len(v) != 4 or v != loaded[0] for v in loaded):
                        raise ValueError('workers disagree on trial sequence')
                    compute_begin = time.perf_counter_ns()
                    futures = [pool.submit(exchange, i, GO, loaded[i], RESULT) for i in range(count)]
                    replies = [f.result() for f in futures]
                    end = time.perf_counter_ns()
                    # Includes QUERY, LOADED, GO and RESULT framing; excludes one-time setup/HELLO.
                    request_bytes = sum(HEADER.size + len(p) + HEADER.size + 4 for p in parts)
                    response_bytes = sum(HEADER.size + 4 + HEADER.size + len(p) for p in replies)
                    prefix = b''.join(U64.pack(v) for v in (upload, end - compute_begin, end - begin,
                                                           request_bytes, response_bytes))
                    send_frame(process.stdin, BATCH_RESULT, prefix + pack_batch(replies))
                elif kind == ROW:
                    row = json.loads(payload)
                    if row.get('correct') != 1 or row.get('mode') != 'tcp_full_protocol':
                        raise ValueError('unverified result')
                    rows.append(row)
                    write_rows(directory / 'raw.csv', rows)
                    print(f'  {identifier(job)} trial {row["trial"]}: Tc={row["client_ms"]:.3f} ms '
                          f'Ts(max)={row["server_max_ms"]:.3f} ms', flush=True)
                elif kind == DONE:
                    if payload or len(rows) != args.trials:
                        raise ValueError('incomplete trial set')
                    break
                else:
                    raise ValueError(f'unknown client frame {kind}')
            process.stdin.close()
            if process.wait(timeout=30):
                raise RuntimeError('C client failed')
            if expired.is_set():
                raise TimeoutError('job exceeded configured deadline')
        status.update(status='complete', samples=len(rows), finished=time.time())
    except BaseException as e:
        status.update(status='interrupted' if isinstance(e, KeyboardInterrupt) else 'failed',
                      error=f'{type(e).__name__}: {e}', samples=len(rows), finished=time.time())
        if expired.is_set():
            status['error'] = 'job exceeded configured deadline'
        raise
    finally:
        if timer:
            timer.cancel()
        if process is not None and process.poll() is None:
            process.kill()
            process.wait()
        for conn in sockets:
            try:
                send_frame(conn, BYE)
            except OSError:
                pass
            conn.close()
        save_json(directory / 'status.json', status)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--profile', choices=['smoke', 'figure2', 'figure3', 'figure4', 'figure5', 'paper'], default='smoke')
    ap.add_argument('--schemes', nargs='+', choices=SCHEMES, default=list(SCHEMES))
    ap.add_argument('--trials', type=int, default=3)
    ap.add_argument('--warmups', type=int, default=1)
    ap.add_argument('--native-repeats', type=int, default=0, help='0: adaptive 100 ms target, capped at 10000 repetitions')
    ap.add_argument('--max-coefficients', type=int, default=5_000_000)
    ap.add_argument('--timeout', type=float, default=86400, help='seconds per job, including HSS key generation')
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()
    if not 1 <= args.trials <= 100000 or not 0 <= args.warmups <= 100 or not 0 <= args.native_repeats <= 100000:
        ap.error('invalid trial/warmup/repetition count')
    if not 1 <= args.max_coefficients <= 2_147_483_647 or args.timeout <= 0:
        ap.error('invalid resource limit')
    config = json.loads(args.config.read_text())
    endpoints = config['servers']
    for i, node in enumerate(endpoints):
        if node.get('server_id') != i or not isinstance(node['host'], str) or not 1 <= node['port'] <= 65535:
            ap.error('servers must be ordered with consecutive server_id values starting at zero')
    if config.get('require_distinct_hosts', True) and len({n['host'] for n in endpoints}) != len(endpoints):
        ap.error('paper deployment requires a distinct host address for every server')
    jobs = jobs_for(args.profile, args.schemes)
    if not jobs:
        ap.error('selection contains no jobs')
    plan = [plan_job(j, args.max_coefficients, endpoints) for j in jobs]
    args.output.mkdir(parents=True, exist_ok=True)
    if args.dry_run:
        save_json(args.output / 'plan.json', plan)
        blocked = [j for j in plan if j['blocked']]
        print(f'{len(plan)} configurations; {len(blocked)} blocked by resource/deployment limits.')
        for j in blocked:
            print(j['id'] + ': ' + '; '.join(j['blocked']))
        return
    if not (ROOT / 'build/vm_client').exists():
        ap.error('run make experiments first')
    verify_build('vm_client')
    source = source_hash()
    settings = {'profile': args.profile, 'schemes': args.schemes, 'trials': args.trials,
                'warmups': args.warmups, 'native_repeats': args.native_repeats,
                'max_coefficients': args.max_coefficients, 'timeout': args.timeout,
                'config': config, 'source_sha256': source}
    fingerprint = hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()
    meta_path = args.output / 'metadata.json'
    if meta_path.exists():
        old = json.loads(meta_path.read_text())
        if not args.resume or old.get('fingerprint') != fingerprint:
            ap.error('output already exists or resume settings/source differ; use a new output directory')
    else:
        save_json(meta_path, {'fingerprint': fingerprint, 'settings': settings, 'client': metadata(),
                             'created': time.time(), 'environment_label': config.get('environment_label', 'unspecified')})
    save_json(args.output / 'plan.json', plan)
    failures = []
    for job, spec in zip(jobs, plan):
        directory = args.output / 'jobs' / spec['id']
        if spec['blocked']:
            save_json(directory / 'status.json', {'status': 'blocked', 'reasons': spec['blocked']})
            failures.append(spec['id'])
            continue
        status_path = directory / 'status.json'
        if args.resume and status_path.exists() and json.loads(status_path.read_text()).get('status') == 'complete':
            continue
        print('Running ' + spec['id'], flush=True)
        try:
            run_job(job, args, config, source, directory)
        except KeyboardInterrupt:
            print('Interrupted; completed jobs are retained. Resume with the same arguments and --resume.', flush=True)
            raise SystemExit(130)
        except Exception as e:
            failures.append(spec['id'])
            print(f'FAILED {spec["id"]}: {e}', flush=True)
    states = {p.parent.name: json.loads(p.read_text())['status'] for p in (args.output / 'jobs').glob('*/status.json')}
    save_json(args.output / 'summary.json', {'jobs': states, 'complete': all(v == 'complete' for v in states.values()),
                                           'finished': time.time()})
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
