#!/usr/bin/env python3
"""Frame-line reduction of the degree-4 mu/nu contraction system (Prop 15.590).

Computes, for any prime p and WITHOUT any Max+/Max- data:
  - the fibers of the frame parameter w (anharmonic S_3 x Frobenius orbits),
  - which fibers are nu-dead (V_4 pairing mechanism, 15.268 generalized),
  - the degree-4 equivariant kernel dimension K_4(p),
  - the pinned functional  sum_f c_f nu_f = V  and hence the RIGOROUS lower
    bound  max_f |nu_hat_f| >= |V| / sum_f |c_f|.

Validated against the four-set implementation at p=5,7,11,13 (kernel dims
1,2,4,6; LB*p^4 = 50.00, 62.36, 91.79, 107.17), ~120x faster; reaches p=31
in under two minutes where the four-set code could not reach p=17 at all.
Note V itself is only defined up to the annihilator normalization; compare
the scale-invariant LB.

Usage:  PS=5,7,11,13,17 python3 scripts/frame_line_system.py

By sharp 3-transitivity of PGL(2,q) every four-set is PGL-equivalent to a frame
{inf,0,1,w}.  For g(z)=(Az+B)/(Cz+D):   g(x)-g(y) = det*(x-y)/((Cx+D)(Cy+D)),
so  C_{gx,gy} = chi(det) * d_x d_y * C_{xy}  with  d_x = chi(Cx+D).
Hence   mu(S) = eps * mu_frame(w),   nu(S) = chi(det) * eps * nu_frame(w),
eps = prod_{x in S} d_x.   Variables: fibers of w under <anharmonic S_3, Frobenius>.
Cost: O(n) per equation instead of enumerating C(n,4) four-sets.

Validates against the four-set implementation (kernel dim, liveDelta, pinned V).
"""
import sys, os, json, time, itertools, random
from fractions import Fraction as F

def make_field(p):
    q=p*p
    r=next(x for x in range(2,p) if pow(x,(p-1)//2,p)==p-1)
    def mul(u,v):
        a1,b1=divmod(u,p); a2,b2=divmod(v,p)
        return p*((a1*a2+r*b1*b2)%p)+((a1*b2+a2*b1)%p)
    def add(u,v):
        a1,b1=divmod(u,p); a2,b2=divmod(v,p)
        return p*((a1+a2)%p)+((b1+b2)%p)
    def neg(u):
        a,b=divmod(u,p); return p*((-a)%p)+((-b)%p)
    def sub(u,v): return add(u,neg(v))
    ONE=p
    inv=[0]*q
    for e in range(1,q):
        x,base,ex=ONE,e,q-2
        while ex:
            if ex&1: x=mul(x,base)
            base=mul(base,base); ex>>=1
        inv[e]=x
    sq=set(mul(x,x) for x in range(1,q))
    chi=[0]*q
    for e in range(1,q): chi[e]= 1 if e in sq else -1
    def frob(e):
        x=e
        for _ in range(p-1): x=mul(x,e)
        return x
    return dict(p=p,q=q,mul=mul,add=add,neg=neg,sub=sub,inv=inv,chi=chi,ONE=ONE,frob=frob)

class Geom:
    def __init__(self,p):
        f=make_field(p)
        self.__dict__.update(f)
        self.INF=self.q               # point at infinity
        self.pts=list(range(self.q))+[self.INF]
    def C(self,x,y):
        if x==y: return 0
        if x==self.INF or y==self.INF: return 1
        return self.chi[self.sub(x,y)]
    # ---- Mobius as (A,B,C,D)
    def apply(self,g,x):
        A,B,Cc,D=g
        if x==self.INF:
            return self.INF if Cc==0 else self.mul(A,self.inv[Cc])
        num=self.add(self.mul(A,x),B); den=self.add(self.mul(Cc,x),D)
        return self.INF if den==0 else self.mul(num,self.inv[den])
    def det(self,g):
        A,B,Cc,D=g
        return self.sub(self.mul(A,D),self.mul(B,Cc))
    def dlift(self,g,x):
        """d_x with d_x d_y C_xy = chi(det) C_{gx,gy}."""
        A,B,Cc,D=g
        dt=self.det(g)
        if Cc==0:
            return self.chi[A] if x==self.INF else self.chi[D]
        pole=self.mul(self.neg(D),self.inv[Cc])
        if x==self.INF: return self.chi[Cc]
        if x==pole:     return self.chi[dt]*self.chi[Cc]
        return self.chi[self.add(self.mul(Cc,x),D)]
    def twist(self,g): return self.chi[self.det(g)]
    def compose(self,g2,g1):
        A2,B2,C2,D2=g2; A1,B1,C1,D1=g1; mul=self.mul; add=self.add
        return (add(mul(A2,A1),mul(B2,C1)), add(mul(A2,B1),mul(B2,D1)),
                add(mul(C2,A1),mul(D2,C1)), add(mul(C2,B1),mul(D2,D1)))
    def minv(self,g):
        A,B,Cc,D=g; return (D,self.neg(B),self.neg(Cc),A)
    def pairing(self,a,b,c,d):
        """Mobius involution with (a b)(c d)."""
        return self.compose(self.minv(self.frame_map(b,a,d)), self.frame_map(a,b,c))
    def frame_map(self,a,b,c):
        """g with g(a)=INF, g(b)=0, g(c)=1."""
        I=self.INF; sub=self.sub; mul=self.mul; neg=self.neg
        if a==I:   return (self.ONE, neg(b), 0, sub(c,b))
        if b==I:   return (0, sub(c,a), self.ONE, neg(a))
        if c==I:   return (self.ONE, neg(b), self.ONE, neg(a))
        A=sub(c,a); B=neg(mul(b,sub(c,a))); Cc=sub(c,b); D=neg(mul(a,sub(c,b)))
        return (A,B,Cc,D)
    def classify(self,S):
        """four-set S (list of 4 points) -> (w, eps, twist) with
           mu(S)=eps*mu_frame(w),  nu(S)=twist*eps*nu_frame(w)."""
        S=sorted(S)
        a,b,c,d=S
        g=self.frame_map(a,b,c)
        w=self.apply(g,d)
        eps=1
        for x in S: eps*=self.dlift(g,x)
        return w,eps,self.twist(g)

def fibers(G):
    """Orbits of L=F_q\\{0,1} under <anharmonic S_3, Frobenius>, with signs.
       Returns lab[w], smu[w], snu[w], nu_dead(set of labels)."""
    q=G.q; ONE=G.ONE; I=G.INF
    L=[w for w in range(q) if w not in (0,ONE)]
    # anharmonic Mobius maps permuting {INF,0,1}
    mobs=[]
    for (a,b,c) in itertools.permutations([I,0,ONE]):
        mobs.append(G.frame_map(a,b,c))
    def step(w):
        out=[]
        for g in mobs:
            w2=G.apply(g,w)
            if w2 in (0,ONE,I): continue
            S=[I,0,ONE,w]
            eps=1
            for x in S: eps*=G.dlift(g,x)
            out.append((w2,eps,G.twist(g)))
        out.append((G.frob(w),1,1))     # Frobenius: preserves C, d=1
        return out
    def pairing_signs(w):
        """V4 stabilizer of the four-SET {INF,0,1,w}: self-maps giving sign
           conditions.  15.268: eps=+1 always; nu dies iff some chi(det)=-1."""
        S=[I,0,ONE,w]; out=[]
        for (x,y) in ((0,1),(0,2),(0,3)):
            rest=[k for k in range(4) if k not in (x,y)]
            tau=G.pairing(S[x],S[y],S[rest[0]],S[rest[1]])
            eps=1
            for z in S: eps*=G.dlift(tau,z)
            out.append((eps,G.twist(tau)))
        return out
    lab={}; smu={}; snu={}; dead=set()
    mu_dead=set()
    for w0 in L:
        if w0 in lab: continue
        lab[w0]=w0; smu[w0]=1; snu[w0]=1
        stack=[w0]
        while stack:
            w=stack.pop()
            for w2,eps,tw in step(w):
                # mu(w)=eps*mu(w2) ; nu(w)=tw*eps*nu(w2)
                nm=smu[w]*eps; nn=snu[w]*tw*eps
                if w2 not in lab:
                    lab[w2]=w0; smu[w2]=nm; snu[w2]=nn; stack.append(w2)
                else:
                    if smu[w2]!=nm: mu_dead.add(w0)
                    if snu[w2]!=nn: dead.add(w0)
    # V4 pairing stabilizer: S -> S, so mu(S)=eps*mu(S), nu(S)=tw*eps*nu(S)
    for w in L:
        for eps,tw in pairing_signs(w):
            if eps==-1: mu_dead.add(lab[w])
            if tw*eps==-1: dead.add(lab[w])
    if mu_dead: raise RuntimeError(f"unexpected mu-dead fibers: {sorted(mu_dead)[:5]}")
    return L,lab,smu,snu,dead

def build(p,extra_pts=14,seed=0,verbose=False):
    t0=time.time()
    G=Geom(p); q=G.q; I=G.INF; ONE=G.ONE; N=1
    L,lab,smu,snu,dead=fibers(G)
    reps=sorted(set(lab.values()))
    live_nu=[r for r in reps if r not in dead]
    iM={r:i for i,r in enumerate(reps)}
    iN={r:len(reps)+i for i,r in enumerate(live_nu)}
    NV=len(reps)+len(live_nu)
    def coefs(S):
        """returns (fiberlabel, mu-sign, nu-sign) for four-set S"""
        w,eps,tw=G.classify(S)
        r=lab[w]
        return r, eps*smu[w], tw*eps*snu[w]
    def addM(vec,S,c):
        r,sm,_=coefs(S); vec[iM[r]]+=c*sm
    def addN(vec,S,c):
        r,_,sn=coefs(S)
        if r not in dead: vec[iN[r]]+=c*sn
    rows=set()
    rng=random.Random(seed)
    pool=[I,0,ONE]+rng.sample([x for x in range(q) if x not in (0,ONE)],
                              min(extra_pts,q-2))
    pool=list(dict.fromkeys(pool))
    # (star mu) and (star nu)
    for i,j,k in itertools.permutations(pool,3):
        if i>j: continue
        cM=[F(0)]*NV; cN=[F(0)]*NV
        for l in G.pts:
            if l in (i,j,k): continue
            c=G.C(k,l)
            if c: addM(cM,[i,j,k,l],c); addN(cN,[i,j,k,l],c)
        rows.add((tuple(cM),F(G.C(i,j))))
        rows.add((tuple(cN),F(-2*G.C(i,k)*G.C(j,k),p)))
    # (out mu) and (out nu): fiber reps x 4 markings, plus random four-sets
    quads=[[I,0,ONE,r] for r in reps]
    for _ in range(3*len(reps)):
        quads.append(rng.sample(G.pts,4))
    for S in quads:
        if len(set(S))<4: continue
        for l_ in S:
            tri=[x for x in S if x!=l_]
            cM=[F(0)]*NV
            for x in G.pts:
                if x in S: continue
                c=G.C(l_,x)
                if c: addM(cM,tri+[x],c)
            addN(cM,S,-p)
            rows.add((tuple(cM),F(0)))
            cN=[F(0)]*NV
            for x in G.pts:
                if x in S: continue
                c=G.C(l_,x)
                if c: addN(cN,tri+[x],c)
            addM(cN,S,-p)
            a_,b_,c_=tri
            corr=(G.C(l_,a_)*G.C(b_,c_)+G.C(l_,b_)*G.C(a_,c_)+G.C(l_,c_)*G.C(a_,b_))
            rows.add((tuple(cN),F(-corr,p)))
    # exact elimination
    M=[list(rw)+[rhs] for rw,rhs in rows]
    rk=0; piv=[]
    for c in range(NV):
        pv=next((i for i in range(rk,len(M)) if M[i][c]!=0),None)
        if pv is None: continue
        M[rk],M[pv]=M[pv],M[rk]; M[rk]=[x/M[rk][c] for x in M[rk]]
        for i2 in range(len(M)):
            if i2!=rk and M[i2][c]!=0: M[i2]=[u-M[i2][c]*v for u,v in zip(M[i2],M[rk])]
        piv.append(c); rk+=1
    incons=any(all(x==0 for x in row[:-1]) and row[-1]!=0 for row in M)
    free=[c for c in range(NV) if c not in piv]
    part=[F(0)]*NV
    for ri,c in enumerate(piv): part[c]=M[ri][NV]
    kvs=[]
    for fc in free:
        kv=[F(0)]*NV; kv[fc]=F(1)
        for ri,c in enumerate(piv): kv[c]=-M[ri][fc]
        kvs.append(kv)
    # pinned functional on the nu block
    off=len(reps); nd=len(live_nu)
    K=[[kv[off+j] for j in range(nd)] for kv in kvs]
    A=[row[:] for row in K]; rk2=0; piv2=[]
    for c in range(nd):
        pv=next((i for i in range(rk2,len(A)) if A[i][c]!=0),None)
        if pv is None: continue
        A[rk2],A[pv]=A[pv],A[rk2]; A[rk2]=[x/A[rk2][c] for x in A[rk2]]
        for i2 in range(len(A)):
            if i2!=rk2 and A[i2][c]!=0: A[i2]=[u-A[i2][c]*v for u,v in zip(A[i2],A[rk2])]
        piv2.append(c); rk2+=1
    free2=[c for c in range(nd) if c not in piv2]
    res=dict(p=p,q=q,n_fibers=len(reps),live_nu=nd,NV=NV,rank=rk,kernel=NV-rk,
             inconsistent=incons,nrows=len(rows),secs=time.time()-t0)
    if len(free2)==1:
        fc=free2[0]
        cvec=[F(0)]*nd; cvec[fc]=F(1)
        for ri,c in enumerate(piv2): cvec[c]=-A[ri][fc]
        V=sum(cvec[j]*part[off+j] for j in range(nd))
        den=sum(abs(x) for x in cvec)
        res.update(V=V,lb=abs(V)/den,lb_p4=float(abs(V)/den)*p**4,ann_dim=1)
    else:
        res.update(ann_dim=len(free2))
    return res

if __name__=="__main__":
    # V is only defined up to scale (annihilator normalization); compare the
    # scale-invariant LB = |V|/sum|c| against the four-set implementation.
    known={5:(1,2,50.00),7:(2,3,62.36),11:(4,5,91.79),13:(6,7,107.17)}
    ps=[int(x) for x in os.environ.get("PS","5,7,11,13").split(",")]
    out=[]
    for p in ps:
        r=build(p)
        tag=""
        if p in known:
            k,ln,lbp=known[p]
            ok=(r['kernel']==k and r['live_nu']==ln
                and abs(r.get('lb_p4',0)-lbp)<0.01)
            tag=f"   vs four-set impl (kernel={k}, liveNu={ln}, LB*p^4={lbp}): {'MATCH' if ok else 'MISMATCH'}"
        print(f"p={p:3d}: fibers={r['n_fibers']:4d} liveNu={r['live_nu']:3d} kernel={r['kernel']:3d} "
              f"V={r.get('V')} LB*p^4={r.get('lb_p4',float('nan')):8.2f} [{r['secs']:.1f}s]{tag}",flush=True)
        out.append({k:(str(v) if isinstance(v,F) else v) for k,v in r.items()})
    json.dump(out,open('/mnt/storage/e1work/leftover3_mu/frame_line.json','w'),indent=1)
