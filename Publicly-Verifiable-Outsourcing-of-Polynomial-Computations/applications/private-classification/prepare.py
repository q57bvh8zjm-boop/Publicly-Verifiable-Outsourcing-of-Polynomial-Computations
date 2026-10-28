"""Rebuild the two fixed public Spambase polynomial classifiers."""
from pathlib import Path
import io, zipfile, json, hashlib, time, sys, platform
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent
SEED = 20261028
Q = 2**252 + 27742317777372353535851937790883648493
SX, SC = 4096, 2**24
threadpool_limits(limits=1)

def metrics(y, score):
    return dict(accuracy=float(accuracy_score(y, score > 0)),
                balanced_accuracy=float(balanced_accuracy_score(y, score > 0)),
                auc=float(roc_auc_score(y, score)))

def timed(fn, repetitions=1000):
    for _ in range(20): fn()
    samples=[]
    for _ in range(7):
        t=time.perf_counter_ns()
        for _ in range(repetitions): fn()
        samples.append((time.perf_counter_ns()-t)/repetitions/1e6)
    return float(np.median(samples))

def quad_export(name, A, l, b, X, y, reference_score):
    m=X.shape[1]; ix=np.rint(X*SX).astype(np.int64)
    # Coefficients of total degree r scaled by SC * SX**(2-r).
    coeff=[int(round(b*SC*SX**2))]+[int(round(v*SC*SX)) for v in l]
    terms=[(-1,-1)]+[(j,-1) for j in range(m)]
    for i in range(m):
        for j in range(i,m):
            coeff.append(int(round((A[i,j] if i==j else A[i,j]+A[j,i])*SC)))
            terms.append((i,j))
    assert max(map(abs,coeff)) < 2**63
    xmax=[int(v) for v in np.max(abs(ix),axis=0)]
    bound=sum(abs(c)*(xmax[i] if i>=0 else 1)*(xmax[j] if j>=0 else 1) for c,(i,j) in zip(coeff,terms))
    assert bound < min(2**126,Q//2)
    # Independent integer oracle for the protocol query subset.
    oracle=[]
    for row in ix[:12]:
        r=[int(v) for v in row]
        oracle.append(sum(c*(r[i] if i>=0 else 1)*(r[j] if j>=0 else 1) for c,(i,j) in zip(coeff,terms)))
    with (ROOT/'data'/f'{name}.fixture').open('w') as f:
        f.write(f'0 {m} 2 {len(coeff)} {len(ix)}\n')
        for c,(i,j) in zip(coeff,terms): f.write(f'{c} {i} {j}\n')
        for label,row in zip(y,ix): f.write(str(int(label))+' '+' '.join(map(str,row))+'\n')
    # Evaluate the compiled quadratic for the native SVM2 timing.
    score=np.einsum('bi,ij,bj->b',X,A,X,optimize=True)+X@l+b
    assert np.max(abs(score-reference_score)) < 1e-6
    baseline_ms=timed(lambda: float(X[0]@A@X[0]+X[0]@l+b))
    np.savez(ROOT/'data'/f'{name}.npz',A=A,l=l,b=b,X=X,y=y,reference_score=score)
    return dict(name=name,kind='quadratic',m=m,d=2,coefficients=len(coeff),
                integer_bound_bits=bound.bit_length(),integer_oracle=[str(v) for v in oracle],
                float_metrics=metrics(y,score),numpy_one_query_ms=baseline_ms,
                quantization={'input_scale':SX,'coefficient_scale':SC},
                model_int64_bytes=8*len(coeff),test_queries=len(X))

def kernel_export(name,model,X,y):
    m=X.shape[1];d=model.degree;sx=1024;sa=2**20
    V=np.rint(model.support_vectors_*sx).astype(np.int64)
    a=np.rint(model.dual_coef_[0]*sa).astype(np.int64)
    ix=np.rint(X*sx).astype(np.int64);offset=m*sx*sx
    bias=int(round(float(model.intercept_[0])*sa*offset**d))
    maxdot=offset+int(np.max(np.sum(abs(V),axis=1)))*int(np.max(abs(ix)))
    bound=abs(bias)+int(np.sum(abs(a)))*maxdot**d
    assert bound < min(2**126,Q//2),(name,bound.bit_length())
    oracle=[]
    for row in ix[:12]:
        dots=V@row+offset
        oracle.append(bias+sum(int(w)*int(z)**d for w,z in zip(a,dots)))
    with (ROOT/'data'/f'{name}.fixture').open('w') as f:
        f.write(f'1 {m} {d} {len(a)} {len(ix)}\n{offset} {bias}\n')
        for w,row in zip(a,V):f.write(str(int(w))+' '+' '.join(map(str,row))+'\n')
        for label,row in zip(y,ix):f.write(str(int(label))+' '+' '.join(map(str,row))+'\n')
    score=model.decision_function(X)
    baseline_ms=timed(lambda: float(model.dual_coef_[0]@((model.support_vectors_@X[0])/m+1)**d+model.intercept_[0]))
    np.savez(ROOT/'data'/f'{name}.npz',V=model.support_vectors_,a=model.dual_coef_[0],bias=model.intercept_[0],X=X,y=y,reference_score=score)
    return dict(name=name,kind='kernel',m=m,d=d,support_vectors=len(a),
                integer_bound_bits=bound.bit_length(),integer_oracle=[str(v) for v in oracle],
                float_metrics=metrics(y,score),numpy_one_query_ms=baseline_ms,
                quantization={'input_and_sv_scale':sx,'alpha_scale':sa},
                model_int64_bytes=8*(len(a)*(m+1))+16,test_queries=len(X))

def svm_experiments(dataset,Xtr,Xte,ytr,yte):
    info=[];m=Xtr.shape[1]
    baselines={}
    for kernel in ['linear','rbf']:
        model=SVC(C=1,kernel=kernel,gamma=1/m).fit(Xtr,ytr)
        baselines[kernel]=metrics(yte,model.decision_function(Xte))
    for d in [2,3]:
        start=time.time()
        model=SVC(C=1,kernel='poly',degree=d,gamma=1/m,coef0=1).fit(Xtr,ytr)
        name=f'{dataset}_svm{d}'
        if d==2:
            V=model.support_vectors_;a=model.dual_coef_[0]
            A=(V.T*a)@V/m**2;l=2*(a@V)/m;b=model.intercept_[0]+a.sum()
            row=quad_export(name,A,l,b,Xte,yte,model.decision_function(Xte))
            row['support_vectors_before_compilation']=len(a)
        else:row=kernel_export(name,model,Xte,yte)
        row.update(training_seconds=time.time()-start,baselines=baselines,
                   train_queries=len(Xtr),task='spam filtering')
        info.append(row);print(json.dumps({k:v for k,v in row.items() if k!='integer_oracle'}),flush=True)
    return info


def main():
    import argparse
    global ROOT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True,help='New output directory')
    parser.add_argument('--dataset',type=Path,default=Path(__file__).resolve().parent/'data/spambase.zip')
    args=parser.parse_args()
    ROOT=args.output.resolve();ROOT.mkdir(parents=True,exist_ok=False)
    (ROOT/'data').mkdir();(ROOT/'results').mkdir()
    with zipfile.ZipFile(args.dataset) as z:
        raw=np.loadtxt(io.BytesIO(z.read('spambase.data')),delimiter=',')
    assert raw.shape==(4601,58) and np.isfinite(raw).all()
    X=np.log1p(raw[:,:-1]);y=raw[:,-1].astype(int)
    tr,te=train_test_split(np.arange(len(y)),test_size=.3,stratify=y,random_state=SEED)
    scaler=StandardScaler().fit(X[tr])
    Xtr=np.clip(scaler.transform(X[tr]),-4,4);Xte=np.clip(scaler.transform(X[te]),-4,4)
    info=svm_experiments('spam',Xtr,Xte,y[tr],y[te])
    result={'seed':SEED,'dataset':'Spambase','split':{'train':tr.tolist(),'test':te.tolist()},
      'workloads':info,'preprocessing':{'operation':'log1p, training-set standardization, clip [-4,4]',
      'mean':scaler.mean_.tolist(),'scale':scaler.scale_.tolist()},
      'dataset_sha256':hashlib.sha256(args.dataset.read_bytes()).hexdigest(),
      'field_order':str(Q),'source':'https://archive.ics.uci.edu/dataset/94/spambase'}
    (ROOT/'results/models.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':main()
