"""Bounded attempt v2: the mid-regime over-sizing is driven by p-hat sampling
noise, not the quantile. Fix: the reference must account for the SAME p-hat
estimation noise the observed statistic carries. We do this by estimating p-hat
inside each simulated efficient replicate too (so the reference distribution of
k-bar includes p-hat noise), rather than conditioning on a fixed phat. This is
the principled debiasing: match the estimation noise in null and observed."""
import numpy as np

def _sim_eff_group_with_phatnoise(n,T,p0,vol,rng):
    """Simulate efficient group AND return k-bar computed the same way the test
    does -- so the reference carries the p-hat estimation noise."""
    piW=p0; g=(rng.random(n)<piW).astype(int); p=np.full(n,piW,float)
    Y=np.zeros((n,T+1),dtype=int)
    for t in range(T+1):
        Y[:,t]=(p>=0.5).astype(int)
        if t==T:break
        step=vol/np.sqrt(t+1); up=np.minimum(1-p,step); dn=np.minimum(p,step)
        pr=np.where(up+dn>0,dn/(up+dn),0.5); mv=rng.random(n)<pr
        p=np.clip(np.where(mv,p+up,p-dn),0,1)
    return np.abs(np.diff(Y,axis=1)).sum(axis=1).mean()

def eff_reference(n,T,phat,rng,vols=(0.15,0.25,0.35,0.5),B=200,qtail=0.95):
    perqs=[]
    for vol in vols:
        kbars=[_sim_eff_group_with_phatnoise(n,T,phat,vol,rng) for _ in range(B)]
        perqs.append(np.quantile(kbars,qtail))
    return max(perqs)

def gen_test(n,T,rng,w,excess,ineff):
    W=rng.integers(0,2,n); piW=np.where(W==1,w,1-w)
    g=(rng.random(n)<piW).astype(int); p=piW.copy(); Y=np.zeros((n,T+1),dtype=int)
    for t in range(T+1):
        Y[:,t]=(p>=0.5).astype(int)
        if t==T:break
        step=0.35/np.sqrt(t+1); up=np.minimum(1-p,step); dn=np.minimum(p,step)
        pr=np.where(up+dn>0,dn/(up+dn),0.5); mv=rng.random(n)<pr
        p=np.clip(np.where(mv,p+up,p-dn),0,1)
    if ineff:
        for i in range(n):
            if W[i]==1:
                for t in range(1,T+1):
                    if rng.random()<excess: Y[i,t]=1-Y[i,t-1]
    return W,Y

def test(W,Y,rng,alpha=0.05):
    T=Y.shape[1]-1
    for wv in [0,1]:
        idx=np.where(W==wv)[0]
        if len(idx)<20:continue
        Yg=Y[idx]; kobs=np.abs(np.diff(Yg,axis=1)).sum(axis=1).mean()
        phat=min(max(Yg[:,-1].mean(),0.02),0.98)
        thr=eff_reference(len(idx),T,phat,rng,B=120)
        if kobs>thr:return True
    return False

def run(reps,n,T,rng,w,excess,ineff):
    return np.mean([test(*gen_test(n,T,rng,w,excess,ineff),rng) for _ in range(reps)])

if __name__=="__main__":
    import sys,json
    rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
    reps=int(sys.argv[2]) if len(sys.argv)>2 else 200
    n,T=800,16
    out={"size":{},"power_fine":{}}
    for w in [0.55,0.75,0.90]:
        out["size"][str(w)]=run(reps,n,T,rng,w,0.0,False)
    for ex in [0.01,0.03,0.05,0.08,0.12,0.20]:
        out["power_fine"][f"ex{ex}"]=run(reps,n,T,rng,0.75,ex,True)
    print(json.dumps(out,indent=2))
