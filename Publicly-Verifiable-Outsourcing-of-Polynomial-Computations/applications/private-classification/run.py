"""Run and verify the public Spambase polynomial-classification experiment."""
from pathlib import Path
import argparse,csv,hashlib,json,platform,statistics,subprocess,time
import numpy as np
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True,help='New measurement directory')
    p.add_argument('--fixtures',type=Path,default=ROOT/'data')
    p.add_argument('--models',type=Path,default=ROOT/'reported/models.json')
    p.add_argument('--binary',type=Path,default=REPO/'build/classification_bench')
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    binary=a.binary.resolve();fixtures=a.fixtures.resolve()
    models={x['name']:x for x in json.loads(a.models.read_text())['workloads']}
    rows=[];quality={};start=time.time()
    for name in ['spam_svm2','spam_svm3']:
      model=models[name]
      for scheme in [1,2]:
        stem=f'{name}-s{scheme}-t1'
        cmd=[str(binary),str(fixtures/(name+'.fixture')),str(scheme),'1','12',str(out/(stem+'-predictions.tsv'))]
        with (out/(stem+'.csv')).open('w') as f,(out/(stem+'.stderr')).open('w') as e:
          subprocess.run(cmd,stdout=f,stderr=e,check=True,timeout=300)
        with (out/(stem+'.csv')).open() as f: raw=list(csv.DictReader(f))
        assert len(raw)==12
        with (out/(stem+'-predictions.tsv')).open() as f:pred=list(csv.DictReader(f,delimiter='\t'))
        assert len(pred)==1381
        scores=[int(x['integer_score']) for x in pred]
        assert list(map(str,scores[:12]))==model['integer_oracle']
        labels=np.array([int(x['label']) for x in pred]);sign=np.array([x>0 for x in scores])
        with np.load(fixtures/(name+'.npz')) as z:
          assert np.array_equal(labels,z['y'])
          changes=int(np.count_nonzero(sign!=(z['reference_score']>0)))
        assert changes==0
        quality[name]={'accuracy':float(np.mean(sign==labels)),'test_queries':len(pred),'float_to_integer_decision_changes':changes}
        row={'name':name,'scheme':raw[0]['scheme'],'t':1,'servers':int(raw[0]['servers']),'trials':len(raw)}
        for k in raw[0]:
          if k.endswith('_ms'):row[k]=statistics.median(float(x[k]) for x in raw)
        for k in ['request_bytes','response_bytes','vk_bytes']:row[k]=int(raw[0][k])
        for k in ['correct','tamper_tests','tamper_rejections']:row[k]=sum(int(x[k]) for x in raw)
        assert row['correct']==12 and row['tamper_tests']==row['tamper_rejections']
        rows.append(row)
        print(name,row['scheme'],'correct',row['correct'],'tamper rejected',row['tamper_rejections'],flush=True)
    total=sum(x['correct'] for x in rows);tamper=sum(x['tamper_tests'] for x in rows)
    assert total==48 and tamper==96
    result={'models':quality,'benchmarks':rows,'correct_protocol_queries':total,'tamper_tests':tamper,'tamper_rejections':tamper}
    (out/'summary.json').write_text(json.dumps(result,indent=2))
    meta={'platform':platform.platform(),'elapsed_seconds':time.time()-start,
      'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
      'source_sha256':{str(f.relative_to(REPO)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [ROOT/'src/classification_bench.c',REPO/'protocols/gamma_fe_sh.c',REPO/'protocols/gamma_ov_sh.c']},
      'fixtures_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in fixtures.glob('*.fixture')},
      'timing':'One warmup, 12 measured queries per configuration; parties run sequentially on one host. No network transport.'}
    (out/'run-metadata.json').write_text(json.dumps(meta,indent=2))
if __name__=='__main__': main()
