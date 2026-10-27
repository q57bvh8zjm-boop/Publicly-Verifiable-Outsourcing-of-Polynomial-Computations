#!/usr/bin/env python3
"""Run measured profiles; save per-trial CSVs, commands, and environment metadata."""
import argparse
import csv
import datetime
import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def output(*args):
    p = subprocess.run(args, cwd=ROOT, text=True, capture_output=True)
    return p.stdout.strip() if p.returncode == 0 else p.stderr.strip()

def source_hashes():
    files = [ROOT / 'Makefile', *ROOT.glob('lib/common/*.[ch]'),
             *ROOT.glob('protocols/*.[ch]'), *ROOT.glob('benchmarks/*.c'),
             *ROOT.glob('scripts/*.py')]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', choices=('smoke', 'paper'), default='smoke')
    parser.add_argument('--trials', type=int, default=3)
    parser.add_argument('--schemes', nargs='+', type=int, choices=(1, 2, 3), default=[1, 2, 3])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.trials < 1:
        parser.error('--trials must be positive')
    destination = args.output.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    subprocess.run(['make', 'benchmark'], cwd=ROOT, check=True)
    runs = []
    if args.profile == 'smoke':
        for scheme in args.schemes:
            runs.extend([(scheme, 2, 2, 1, 'dense'), (scheme, 2, 2, 1, 'var')])
        for scheme in args.schemes:
            if scheme != 3:
                runs.extend((scheme, 2, d, 1, 'dense') for d in (4, 8, 10))
    else:
        for scheme in args.schemes:
            runs.extend((scheme, m, 2, 1, 'dense') for m in range(200, 2001, 200))
            runs.extend((scheme, m, d, 1, 'dense') for d in (4, 8, 10) for m in range(1, 16))
            if scheme != 3:
                runs.extend((scheme, m, 2, t, 'dense') for t in (2, 3) for m in range(200, 2001, 200))
    metadata = {'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'platform': platform.platform(), 'machine': platform.machine(),
                'processor': platform.processor(), 'cpu_count': os.cpu_count(),
                'compiler': output(os.environ.get('CC', 'cc'), '--version'),
                'revision': output('git', 'rev-parse', 'HEAD'),
                'working_tree': output('git', 'status', '--short'),
                'source_sha256': source_hashes(), 'profile': args.profile,
                'trials': args.trials, 'commands': [], 'status': 'running'}
    meta_path = destination / 'metadata.json'
    all_rows = []
    try:
        for index, (scheme, m, d, t, workload) in enumerate(runs):
            command = ['./build/benchmark', str(scheme), str(m), str(d), str(t), str(args.trials), workload]
            metadata['commands'].append(command)
            meta_path.write_text(json.dumps(metadata, indent=2) + '\n')
            print('Running:', ' '.join(command), flush=True)
            with (destination / f'{index:03d}.csv').open('w') as out, (destination / f'{index:03d}.stderr').open('w') as err:
                subprocess.run(command, cwd=ROOT, stdout=out, stderr=err, check=True)
            with (destination / f'{index:03d}.csv').open() as source:
                rows = list(csv.DictReader(source))
                if len(rows) != args.trials or any(r['correct'] != '1' for r in rows):
                    raise RuntimeError('Incomplete or incorrect trial output')
                all_rows.extend(rows)
        with (destination / 'raw.csv').open('w', newline='') as out:
            writer = csv.DictWriter(out, fieldnames=list(all_rows[0]))
            writer.writeheader(); writer.writerows(all_rows)
        metadata['status'] = 'complete'
    finally:
        metadata['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        meta_path.write_text(json.dumps(metadata, indent=2) + '\n')

if __name__ == '__main__':
    main()
