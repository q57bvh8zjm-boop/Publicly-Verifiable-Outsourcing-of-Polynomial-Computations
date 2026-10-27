#!/usr/bin/env python3
"""One single-threaded C computation worker per independently deployed VM."""
import argparse
import json
import signal
import socket
import subprocess
import threading
import time
from vm_protocol import ROOT, VERSION, HELLO, HELLO_REPLY, ERROR, metadata, recv_frame, send_frame, verify_build


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--server-id', type=int, required=True)
    ap.add_argument('--bind', default='127.0.0.1')
    ap.add_argument('--port', type=int, default=9100)
    ap.add_argument('--max-coefficients', type=int, default=5_000_000)
    ap.add_argument('--job-timeout', type=float, default=86400)
    args = ap.parse_args()
    if not 0 <= args.server_id < 128 or not 1 <= args.max_coefficients <= 2_147_483_647 or args.job_timeout <= 0:
        ap.error('invalid server id, coefficient limit, or timeout')
    worker = ROOT / 'build/vm_worker'
    if not worker.exists():
        ap.error('run make experiments first')
    lock = threading.Lock()
    children, children_lock = set(), threading.Lock()

    def handle(conn, address):
        acquired, process = False, None
        try:
            conn.settimeout(30)
            _, payload = recv_frame(conn, HELLO)
            request = json.loads(payload)
            if request.get('wire_version') != VERSION or request.get('server_id') != args.server_id:
                raise ValueError('wire version or server identity mismatch')
            acquired = lock.acquire(timeout=5)
            if not acquired:
                raise RuntimeError('this VM already has an active experiment')
            verify_build('vm_worker')
            info = metadata()
            info.update(server_id=args.server_id, max_coefficients=args.max_coefficients,
                        job_timeout=args.job_timeout, connected_at=time.time())
            send_frame(conn, HELLO_REPLY, json.dumps(info).encode())
            # Restore blocking mode before the socket is inherited as C stdin/stdout.
            conn.settimeout(None)
            conn.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            process = subprocess.Popen([str(worker), str(args.server_id), str(args.max_coefficients)],
                                       stdin=conn, stdout=conn)
            with children_lock:
                children.add(process)
            process.wait(timeout=args.job_timeout)
            if process.returncode:
                print(f'{address}: worker failed with exit code {process.returncode}', flush=True)
        except Exception as e:
            print(f'{address}: {type(e).__name__}: {e}', flush=True)
            try:
                send_frame(conn, ERROR, str(e).encode())
            except OSError:
                pass
        finally:
            if process is not None and process.poll() is None:
                process.kill()
                process.wait()
            with children_lock:
                children.discard(process)
            conn.close()
            if acquired:
                lock.release()

    def stop(_signum, _frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, stop)
    try:
        with socket.socket() as listener:
            listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            listener.bind((args.bind, args.port))
            listener.listen(16)
            print(f'server {args.server_id} listening on {args.bind}:{args.port}; one active C worker', flush=True)
            while True:
                conn, address = listener.accept()
                threading.Thread(target=handle, args=(conn, address), daemon=True).start()
    except KeyboardInterrupt:
        pass
    finally:
        with children_lock:
            active = list(children)
        for process in active:
            if process.poll() is None:
                process.kill()
        for process in active:
            process.wait()


if __name__ == '__main__':
    main()
