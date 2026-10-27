#!/usr/bin/env python3
"""Measure Table 4 primitives on the actual client VM (100 fresh trials by default)."""
import argparse
import csv
from pathlib import Path
import statistics
import subprocess
from vm_protocol import ROOT, metadata, save_json, verify_build


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--trials', type=int, default=100)
    ap.add_argument('--inner-repeats', type=int, default=1000)
    ap.add_argument('--skip-paillier', action='store_true', help='quick validation only; leaves Table 4 incomplete')
    args = ap.parse_args()
    if not 1 <= args.trials <= 100000 or not 1 <= args.inner_repeats <= 1000000:
        ap.error('invalid repetition count')
    if args.output.exists():
        ap.error('output exists; use a fresh directory')
    verify_build('microbenchmark')
    args.output.mkdir(parents=True)
    command = [str(ROOT / 'build/microbenchmark'), str(args.trials), str(args.inner_repeats), str(int(not args.skip_paillier))]
    save_json(args.output / 'metadata.json', {'client': metadata(), 'command': command,
              'group_operation': 'random-base Ristretto255 scalar multiplication',
              'paillier': 'one randomized encryption under a fresh >=3072-bit HSS-compatible key; not complete HSS Input'})
    with (args.output / 'raw.csv').open('w') as out, (args.output / 'stderr.log').open('w') as err:
        result = subprocess.run(command, stdout=out, stderr=err)
    if result.returncode:
        save_json(args.output / 'status.json', {'status': 'failed', 'returncode': result.returncode})
        raise SystemExit(result.returncode)
    groups = {}
    with (args.output / 'raw.csv').open() as f:
        for row in csv.DictReader(f):
            groups.setdefault(row['operation'], []).append(float(row['ns_per_operation']))
    save_json(args.output / 'summary.json', {k: {'samples': len(v), 'mean_ns': statistics.mean(v),
              'median_ns': statistics.median(v), 'stdev_ns': statistics.stdev(v) if len(v) > 1 else 0}
              for k, v in groups.items()})
    save_json(args.output / 'status.json', {'status': 'complete', 'table_complete': not args.skip_paillier})


if __name__ == '__main__':
    main()
