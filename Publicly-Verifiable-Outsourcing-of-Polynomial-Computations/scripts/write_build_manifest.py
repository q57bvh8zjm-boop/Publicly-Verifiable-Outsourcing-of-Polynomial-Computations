#!/usr/bin/env python3
"""Record the source, build options, and binaries produced by make experiments."""
import argparse
import hashlib
import os
from vm_protocol import ROOT, save_json, source_hash

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--compiler', required=True)
ap.add_argument('--cflags', required=True)
ap.add_argument('--libs', required=True)
args = ap.parse_args()
save_json(ROOT / 'build/experiment-build.json', {
    'source_sha256': source_hash(), 'compiler': args.compiler, 'cflags': args.cflags,
    'libs': args.libs, 'dependency_prefix': os.environ.get('PV_PREFIX'),
    'binaries': {name: hashlib.sha256((ROOT / 'build' / name).read_bytes()).hexdigest()
                 for name in ('vm_client', 'vm_worker', 'microbenchmark')},
})
