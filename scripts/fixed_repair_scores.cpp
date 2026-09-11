// Exact projective Gray-code scoring of one prescribed repair/addition lattice.
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <omp.h>
int main(int argc, char** argv) {
  if (argc != 2) return 2;
  int workers=std::atoi(argv[1]);
  if(workers<1 || workers>omp_get_num_procs()) return 2;
  int k,N,r; if(!(std::cin>>k>>N>>r)) return 2;
  if(k<2 || N<k || N>24 || r<0 || N-k+r>16) return 2;
  std::vector<int> B(N*N), u(r),v(r);
  for(auto &a:B) if(!(std::cin>>a)) return 2;
  for(int j=0;j<r;j++) if(!(std::cin>>u[j]>>v[j]) || u[j]<0 || u[j]>=v[j] || v[j]>=k) return 2;
  for(int i=0;i<N;i++) for(int j=0;j<N;j++)
    if(B[i*N+j]!=B[j*N+i] || (i==j ? B[i*N+j]!=0 : std::abs(B[i*N+j])!=1)) return 2;
  int m=N-k, total=1<<(m+r);
  std::vector<int> lo(total),hi(total),wm(total),ns(total);
  std::vector<long long> sums(total),squares(total);
  #pragma omp parallel for num_threads(workers) schedule(dynamic)
  for(int state=0;state<total;state++) {
    std::vector<int> ids;
    for(int i=0;i<k;i++) ids.push_back(i);
    for(int i=0;i<m;i++) if(state&(1<<i)) ids.push_back(k+i);
    int n=ids.size(); ns[state]=n;
    int A[24][24]={}, x[24],h[24]={};
    for(int i=0;i<n;i++) { x[i]=1; for(int j=0;j<n;j++) A[i][j]=B[ids[i]*N+ids[j]]; }
    for(int j=0;j<r;j++) if(!(state&(1<<(m+j)))) A[u[j]][v[j]]=A[v[j]][u[j]]=-A[u[j]][v[j]];
    int q=0;
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) {h[i]+=A[i][j]; if(i<j) q+=A[i][j];}
    int mn=q,mx=q,w=0,mask=0;
    lo[state]=q; hi[state]=q;
    long long sum=q,sq=1LL*q*q;
    for(int t=1;t<(1<<(n-1));t++) {
      int j=__builtin_ctz((unsigned)t)+1,old=x[j];
      q-=2*old*h[j]; x[j]=-old; mask^=1<<(j-1);
      for(int i=0;i<n;i++) h[i]-=2*old*A[i][j];
      mn=std::min(mn,q); mx=std::max(mx,q);
      if(std::abs(q)>std::max(std::abs(lo[state]),std::abs(hi[state]))) {wm[state]=mask;}
      lo[state]=mn;hi[state]=mx;
      sum+=q; sq+=1LL*q*q;
    }
    // Record a witness via a second, independent direct check of the extrema mask.
    // The initial all-plus state may itself attain the maximum.
    w=wm[state]; int direct=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) {
      int xi=i==0?1:((w>>(i-1)&1)?-1:1),xj=(w>>(j-1)&1)?-1:1;
      direct+=A[i][j]*xi*xj;
    }
    if(std::abs(direct)!=std::max(-mn,mx)) std::abort();
    lo[state]=mn; hi[state]=mx;sums[state]=sum;squares[state]=sq;
  }
  for(int s=0;s<total;s++) {
    long long spins=1LL<<(ns[s]-1);
    if(sums[s]!=0 || squares[s]!=spins*ns[s]*(ns[s]-1)/2) return 3;
    std::cout<<s<<' '<<ns[s]<<' '<<lo[s]<<' '<<hi[s]<<' '<<wm[s]<<'\n';
  }
}
