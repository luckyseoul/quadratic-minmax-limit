// defect: decide whether the Paley conference matrix C_{q+1} (q prime, q = 1 mod 4) has a Boolean x whose
// local fields h_i = x_i (Cx)_i form a prescribed multiset  { c : n - s } u { c + 4 t : count }  (t != 0).
//
// Reduction (proved in the note).  Put e = (Cx - c x)/4, i.e. e_i = t_i x_i.  Then (C + c) e = m x with
// m = (q - c^2)/4, and conversely any integer e with (C + c) e = m x for a Boolean x, e_i = t_i x_i, gives
// Cx = c x + 4 e.  So x exists iff there are a support S, values t_i on S and signs x_i on S with
//      (Ce)_i = x_i (m - c t_i)  for i in S,          |(Ce)_i| = m  for i not in S.
// Symmetry: PSL(2,q) acts 2-transitively on the points by signed permutations preserving C, and x -> -x is a
// symmetry, so the first support vertex is pinned to the point inf with x = +1 and the second to the point 0.
// The stabiliser of the ordered pair (inf, 0) is the split torus j -> s j (s a square), with two orbits on the
// remaining points, so the third support vertex is pinned to the point 1 or to the least nonsquare nu.
// Later vertices of one role (one t value) are added in increasing index order.
//
// usage: defect q c  t1 n1  t2 n2 ...   (roles in the order they are to be placed; all solutions are counted)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int q,n,c,m; static signed char *C;
static int nr, rt[16], rn[16];          // roles: t value, count
static int S[256], Sx[256], St[256], ns; // support so far
static int inS[256]; static long p[256]; // p = C e (partial)
static int total_w, placed_w; static unsigned long long nodes=0, sols=0; static int first_printed=0;
static int last_idx_role[16], cnt_role[16];
static int seq_role[256], nseq; static int nu, second_sign=0;

static int feasible(int rem){ // rem = remaining weight still to be placed
  for(int i=0;i<n;i++){
    long pi=p[i];
    if(inS[i]>=0){ int k=inS[i]; long tgt=(long)Sx[k]*(m-(long)c*St[k]); long d=tgt-pi; if(d<0)d=-d;
      if(d>rem || ((d-rem)&1)) return 0; }
    else{
      int ok=0; long d1=pi-m; if(d1<0)d1=-d1; long d2=pi+m; if(d2<0)d2=-d2;
      if((d1<=rem && !((d1-rem)&1)) || (d2<=rem && !((d2-rem)&1))) ok=1;
      if(!ok){ // could i still join the support in a role with vertices left?
        for(int r=0;r<nr && !ok;r++) if(cnt_role[r]<rn[r]){ int at=rt[r]<0?-rt[r]:rt[r]; int rr=rem-at; if(rr<0) continue;
          for(int sx=-1;sx<=1 && !ok;sx+=2){ long tgt=(long)sx*(m-(long)c*rt[r]); long d=tgt-pi; if(d<0)d=-d; if(d<=rr && !((d-rr)&1)) ok=1; } }
      }
      if(!ok) return 0;
    }
  }
  return 1;
}
static void add(int j,int x,int t){ long e=(long)t*x; const signed char*col=C+(size_t)j*n; for(int i=0;i<n;i++) p[i]+=col[i]*e; inS[j]=ns; S[ns]=j; Sx[ns]=x; St[ns]=t; ns++; }
static void del(void){ ns--; int j=S[ns]; long e=(long)St[ns]*Sx[ns]; const signed char*col=C+(size_t)j*n; for(int i=0;i<n;i++) p[i]-=col[i]*e; inS[j]=-1; }
static void dfs(int k){ // place the k-th support vertex
  nodes++;
  int rem=total_w-placed_w;
  if(k==nseq){ // complete: all rows already exact (rem=0 feasibility)
    sols++;
    if(!first_printed){ first_printed=1; printf("SOLUTION x="); for(int i=0;i<n;i++){ int xi= inS[i]>=0? Sx[inS[i]] : (p[i]>0?1:-1); putchar(xi>0?'+':'-'); } printf("\n"); fflush(stdout); if(getenv("FIRST")){ printf("q=%d c=%d first solution after %llu nodes\n",q,c,nodes); exit(0);} } 
    return; }
  int r=seq_role[k], t=rt[r], at=t<0?-t:t;
  int lo, hi; int cand3[2];
  if(k==0){ lo=0; hi=0; } else if(k==1){ lo=1; hi=1; }
  else if(k==2){ cand3[0]=2; cand3[1]=1+nu; lo=0; hi=1; }   // third vertex: the point 1 or the nonsquare nu
  else { lo=last_idx_role[r]+1; if(lo<2) lo=2; hi=n-1-(rn[r]-cnt_role[r]-1); }
  for(int jj=lo;jj<=hi;jj++){
    int j = (k==2)? cand3[jj] : jj;
    if(inS[j]>=0) continue;
    for(int x=1;x>=-1;x-=2){
      if(k==0 && x<0) continue;
      if(k==1 && second_sign && x!=second_sign) continue;
      // quick local test: row j itself must be able to reach its target
      { long tgt=(long)x*(m-(long)c*t); long d=tgt-p[j]; if(d<0)d=-d; int rr=rem-at; if(d>rr || ((d-rr)&1)) continue; }
      int save=last_idx_role[r];
      add(j,x,t); placed_w+=at; cnt_role[r]++; if(k>=3) last_idx_role[r]=j;
      if(feasible(total_w-placed_w)) dfs(k+1);
      last_idx_role[r]=save; cnt_role[r]--; placed_w-=at; del();
    }
  }
}
int main(int argc,char**argv){
  q=atoi(argv[1]); c=atoi(argv[2]); n=q+1; m=(q-c*c)/4; if((q-c*c)%4){fprintf(stderr,"bad c\n");return 1;}
  nr=(argc-3)/2; total_w=0; nseq=0;
  for(int r=0;r<nr;r++){ rt[r]=atoi(argv[3+2*r]); rn[r]=atoi(argv[4+2*r]); int at=rt[r]<0?-rt[r]:rt[r]; total_w+=at*rn[r]; for(int k=0;k<rn[r];k++) seq_role[nseq++]=r; last_idx_role[r]=-1; cnt_role[r]=0; }
  C=malloc((size_t)n*n); int *chi=calloc(q,sizeof(int)); for(int a=1;a<q;a++) chi[(int)(((long long)a*a)%q)]=1; for(int a=1;a<q;a++) if(!chi[a]) chi[a]=-1;
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){int v; if(i==j)v=0; else if(i==0||j==0)v=1; else v=chi[((i-j)%q+q)%q]; C[(size_t)i*n+j]=v;}
  for(int i=0;i<n;i++){inS[i]=-1;p[i]=0;} ns=0; placed_w=0;
  nu=2; while(chi[nu]!=-1) nu++;
  if(getenv("SECOND_SIGN")) second_sign=atoi(getenv("SECOND_SIGN"));
  dfs(0);
  printf("q=%d c=%d m=%d support=%d weight=%d nodes=%llu solutions_found=%llu\n",q,c,m,nseq,total_w,nodes,sols);
  return 0;
}
