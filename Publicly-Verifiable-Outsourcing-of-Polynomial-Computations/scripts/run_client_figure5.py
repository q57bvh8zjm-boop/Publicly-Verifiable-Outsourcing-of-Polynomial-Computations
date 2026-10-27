#!/usr/bin/env python3
"""Reproduce Figure 5 from measured client calls and valid reference responses."""
import argparse, csv, datetime, hashlib, json, os, platform, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def capture(*args):
    p=subprocess.run(args,cwd=ROOT,text=True,capture_output=True)
    return p.stdout.strip() if p.returncode==0 else p.stderr.strip()
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--trials',type=int,default=3)
    args=parser.parse_args()
    if args.trials<1: parser.error('--trials must be positive')
    out=args.output.resolve(); out.mkdir(parents=True,exist_ok=False)
    subprocess.run(['make','client-benchmark'],cwd=ROOT,check=True)
    source=[ROOT/'Makefile',*ROOT.glob('lib/common/*.[ch]'),*ROOT.glob('protocols/*.[ch]'),*ROOT.glob('benchmarks/*.c'),*ROOT.glob('scripts/*.py')]
    metadata={'started_utc':now(),'status':'running','platform':platform.platform(),'machine':platform.machine(),
      'cpu_model':capture('/usr/sbin/sysctl','-n','machdep.cpu.brand_string'),'memory_bytes':capture('/usr/sbin/sysctl','-n','hw.memsize'),
      'cpu_count':os.cpu_count(),'compiler':capture('cc','--version'),'revision':capture('git','rev-parse','HEAD'),
      'working_tree':capture('git','status','--short'),
      'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source)},
      'trials':args.trials,'warmup_trials_per_configuration':1,'verify_decode_inner_runs':25,
      'native_timing':'one warmup Eval, then 1..10000 iterations targeting 20 ms; includes result cleanup',
      'mode':'client-only microbenchmark; valid reference responses; no server execution or network timing',
      'hss_setup':'one fresh key per degree, provisioned for m=15; setup excluded from client time',
      'aggregation':'median per-trial client/native ratio; log10 applied after median',
      'commands':[],'completed_degrees':[]}
    meta=out/'metadata.json'
    def save(): meta.write_text(json.dumps(metadata,indent=2)+'\n')
    save(); rows=[]
    try:
        for d in (4,8,10):
            command=['./build/client_benchmark',str(d),'15',str(args.trials),'1']; metadata['commands'].append(command);save()
            print(now(),'Running',' '.join(command),flush=True)
            with (out/f'd{d}.csv').open('w') as stdout,(out/f'd{d}.stderr').open('w') as stderr:
                subprocess.run(command,cwd=ROOT,stdout=stdout,stderr=stderr,check=True)
            with (out/f'd{d}.csv').open() as f: current=list(csv.DictReader(f))
            assert len(current)==45*args.trials and all(r['correct']=='1' for r in current)
            rows+=current; metadata['completed_degrees'].append(d);save()
        with (out/'raw.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        metadata['status']='complete'
    except BaseException:
        metadata['status']='failed';raise
    finally:
        metadata['finished_utc']=now();save()
if __name__=='__main__':main()
