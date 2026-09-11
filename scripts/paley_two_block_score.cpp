// Exact scorer for K(B)=[[C,B],[B^T,-C]], where C is q=13 Paley.
// B is read as 196 whitespace-separated signs.  The global cube symmetry
// fixes the first left spin, leaving 2^27 states, all evaluated exactly.
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <vector>

struct Witness { int value; int energy; int xmask; int ymask; };

int main(int argc, char** argv) {
  if (argc != 2) { std::cerr << "usage: scorer B.txt\n"; return 2; }
  constexpr int q=13, n=14, xs=1<<(n-1), ys=1<<n;
  int chi[q]={};
  for (int a=1;a<q;++a) chi[(a*a)%q]=1;
  for (int a=1;a<q;++a) if (!chi[a]) chi[a]=-1;
  int C[n][n]={};
  for (int j=1;j<n;++j) C[0][j]=C[j][0]=1;
  for (int i=0;i<q;++i) for (int j=i+1;j<q;++j)
    C[i+1][j+1]=C[j+1][i+1]=chi[(j-i)%q];
  int B[n][n]; std::ifstream in(argv[1]);
  for (int i=0;i<n;++i) for (int j=0;j<n;++j)
    if (!(in>>B[i][j]) || (B[i][j]!=1 && B[i][j]!=-1)) { std::cerr<<"invalid B\n"; return 2; }
  std::vector<int8_t> x(xs*n), y(ys*n);
  std::vector<int16_t> qx(xs), qy(ys);
  int source=0;
  for(int u=0;u<xs;++u) { auto z=&x[u*n]; z[0]=1;
    for(int i=1;i<n;++i) z[i]=((u>>(i-1))&1)?1:-1;
    int e=0; for(int i=0;i<n;++i) for(int j=i+1;j<n;++j)e+=C[i][j]*z[i]*z[j];
    qx[u]=e; source=std::max(source,std::abs(e)); }
  for(int u=0;u<ys;++u) { auto z=&y[u*n];
    for(int i=0;i<n;++i) z[i]=((u>>i)&1)?1:-1;
    int e=0; for(int i=0;i<n;++i) for(int j=i+1;j<n;++j)e+=C[i][j]*z[i]*z[j]; qy[u]=e; }
  std::atomic<int> best(0); Witness witness{0,0,0,0}; std::mutex wlock;
  #pragma omp parallel
  {
    int local_best=0, local_energy=0, lx=0, ly=0;
    std::vector<int16_t> h(n);
    #pragma omp for schedule(static)
    for(int v=0;v<ys;++v) {
      for(int i=0;i<n;++i) { int s=0; for(int j=0;j<n;++j)s+=B[i][j]*y[v*n+j]; h[i]=s; }
      for(int u=0;u<xs;++u) {
        int cross=0; for(int i=0;i<n;++i) cross+=x[u*n+i]*h[i];
        int value=std::abs((int)qx[u]-(int)qy[v]+cross);
        if(value>local_best) {local_best=value;local_energy=(int)qx[u]-(int)qy[v]+cross;lx=u;ly=v;}
      }
    }
    int old=best.load(std::memory_order_relaxed);
    while(local_best>old && !best.compare_exchange_weak(old,local_best,std::memory_order_relaxed)) {}
    if(local_best>=best.load(std::memory_order_relaxed)) { std::lock_guard<std::mutex> lock(wlock);
      if(local_best>witness.value) witness={local_best,local_energy,lx,ly}; }
  }
  // A second exact pass retains a small active band.  It is used only for
  // ranking integral repair proposals; the first pass already establishes
  // the reported Boolean norm.
  const int cutoff = best.load() - 2;
  std::vector<Witness> active;
  std::mutex alock;
  #pragma omp parallel
  {
    std::vector<Witness> local;
    std::vector<int16_t> h(n);
    #pragma omp for schedule(static)
    for (int v=0; v<ys; ++v) {
      for (int i=0;i<n;++i) { int s=0; for (int j=0;j<n;++j)s+=B[i][j]*y[v*n+j]; h[i]=s; }
      for (int u=0;u<xs;++u) {
        int cross=0; for (int i=0;i<n;++i) cross+=x[u*n+i]*h[i];
        int energy=(int)qx[u]-(int)qy[v]+cross;
        if (std::abs(energy)>cutoff && local.size()<512)
          local.push_back({std::abs(energy),energy,u,v});
      }
    }
    std::lock_guard<std::mutex> lock(alock);
    active.insert(active.end(), local.begin(), local.end());
  }
  std::sort(active.begin(),active.end(),[](const Witness&a,const Witness&b){return a.value>b.value;});
  if (active.size()>512) active.resize(512);
  std::cout << "{\"n\":28,\"source_phi\":"<<source<<",\"phi\":"<<best.load()
            <<",\"energy\":"<<witness.energy<<",\"xmask\":"<<witness.xmask<<",\"ymask\":"<<witness.ymask
            <<",\"workers\":"<<omp_get_max_threads()<<",\"active\":[";
  for (size_t k=0;k<active.size();++k) { if(k) std::cout<<','; const auto&a=active[k];
    std::cout<<"{\"value\":"<<a.value<<",\"energy\":"<<a.energy<<",\"xmask\":"<<a.xmask<<",\"ymask\":"<<a.ymask<<'}'; }
  std::cout << "]}\n";
}
