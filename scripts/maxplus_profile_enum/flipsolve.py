"""Exact flip-assignment solver: DFS with per-point interval feasibility,
incremental line-local updates, and weight bookkeeping. No 2^n brute."""
import numpy as np

def flip_solve(p,k,sigstd,idxs,fvec,Tm,thi,tlo,out,node_cap=8_000_000):
    q=Tm.shape[1]
    A=np.zeros(q,dtype=np.int64)
    for j in range(k): A=A+sigstd[j][Tm[j]]
    diff=A-thi
    if (diff%(2*p)).any(): return
    g=diff//(2*p)                      # Xi(x) must be in {g, g+1}
    if (g>k).any() or (g+1<0).any(): return
    bits=[(j,int(s)) for j in range(k) for s in idxs[j]]
    nb=len(bits)
    # per-point: assigned Xi so far, and #free bits covering the point
    assigned=np.zeros(q,dtype=np.int64)
    freec=np.zeros(q,dtype=np.int64)
    lines={}
    for j,s in bits:
        line=np.where(Tm[j]==s)[0]
        lines[(j,s)]=line
        freec[line]+=1
    # quick global feasibility
    if ((assigned>g+1)|(assigned+freec<g)).any(): return
    remw=[int(f) for f in fvec]        # remaining weight per profile
    remn=[len(idxs[j]) for j in range(k)]
    nodes=[0]
    sols=[]
    def dfs(i):
        if nodes[0]>node_cap: raise RuntimeError("flip node cap")
        nodes[0]+=1
        if i==nb:
            Xi=assigned
            y=np.where(Xi==g,1,-1).astype(np.int8)
            # sanity: Xi in {g,g+1}
            sols.append(y.copy())
            return
        j,s=bits[i]
        line=lines[(j,s)]
        for val in (0,1):
            if val==1 and remw[j]==0: continue
            if val==0 and remn[j]-1<remw[j]: continue
            assigned[line]+=val; freec[line]-=1
            remw[j]-=val; remn[j]-=1
            bad=((assigned[line]>g[line]+1)|(assigned[line]+freec[line]<g[line])).any()
            if not bad: dfs(i+1)
            assigned[line]-=val; freec[line]+=1
            remw[j]+=val; remn[j]+=1
    dfs(0)
    out.extend(sols)
