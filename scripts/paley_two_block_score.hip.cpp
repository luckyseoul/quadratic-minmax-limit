// HIP exact full-cube scorer for K(B)=[[C,B],[B^T,-C]], q=13 Paley C.
// This is an accelerator parity backend; it evaluates all 2^27 projective
// states using integer arithmetic and reports one exact maximizing witness.
#include <hip/hip_runtime.h>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>

constexpr int N=14, XS=1<<13, YS=1<<14;
constexpr unsigned long long TOTAL=(unsigned long long)XS*YS;
__constant__ int8_t DC[N*N];
__constant__ int8_t DB[N*N];
#define HIP_CHECK(call) do { hipError_t e=(call); if(e!=hipSuccess){std::cerr<<hipGetErrorString(e)<<"\n";return 1;} } while(0)

__device__ inline int sx(unsigned mask,int i) { return i==0 ? 1 : ((mask>>(i-1))&1 ? 1 : -1); }
__device__ inline int sy(unsigned mask,int i) { return (mask>>i)&1 ? 1 : -1; }
__device__ inline unsigned long long packed(unsigned maskx,unsigned masky) {
  int ex=0, ey=0, cross=0;
  for(int i=0;i<N;++i) for(int j=i+1;j<N;++j) {
    ex += DC[i*N+j]*sx(maskx,i)*sx(maskx,j);
    ey += DC[i*N+j]*sy(masky,i)*sy(masky,j);
  }
  for(int i=0;i<N;++i) for(int j=0;j<N;++j) cross += DB[i*N+j]*sx(maskx,i)*sy(masky,j);
  int energy=ex-ey+cross, value=energy<0?-energy:energy;
  unsigned long long index=((unsigned long long)masky<<13)|maskx;
  return ((unsigned long long)(unsigned)value<<32)|index;
}
__global__ void score(unsigned long long* answer) {
  __shared__ unsigned long long blockbest;
  if(threadIdx.x==0) blockbest=0;
  __syncthreads();
  unsigned long long local=0;
  for(unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x;
      i<TOTAL; i+=(unsigned long long)gridDim.x*blockDim.x) {
    unsigned long long v=packed((unsigned)(i&(XS-1)),(unsigned)(i>>13));
    if(v>local) local=v;
  }
  atomicMax(&blockbest,local);
  __syncthreads();
  if(threadIdx.x==0) atomicMax(answer,blockbest);
}
int host_energy(const int8_t* C,const int8_t* B,unsigned x,unsigned y) {
  int ex=0,ey=0,cross=0;
  auto X=[x](int i){return i==0?1:((x>>(i-1))&1?1:-1);};
  auto Y=[y](int i){return (y>>i)&1?1:-1;};
  for(int i=0;i<N;++i)for(int j=i+1;j<N;++j){ex+=C[i*N+j]*X(i)*X(j);ey+=C[i*N+j]*Y(i)*Y(j);}
  for(int i=0;i<N;++i)for(int j=0;j<N;++j)cross+=B[i*N+j]*X(i)*Y(j);
  return ex-ey+cross;
}
int main(int argc,char**argv) {
  if(argc!=2){std::cerr<<"usage: hip_scorer B.txt\n";return 2;}
  int8_t C[N*N]={},B[N*N],chi[13]={};
  for(int a=1;a<13;++a)chi[(a*a)%13]=1;
  for(int a=1;a<13;++a)if(!chi[a])chi[a]=-1;
  for(int j=1;j<N;++j)C[j]=C[j*N]=1;
  for(int i=0;i<13;++i)for(int j=i+1;j<13;++j)C[(i+1)*N+j+1]=C[(j+1)*N+i+1]=chi[(j-i)%13];
  std::ifstream in(argv[1]); int z;
  for(auto &v:B){if(!(in>>z)||(z!=1&&z!=-1)){std::cerr<<"invalid B\n";return 2;}v=(int8_t)z;}
  HIP_CHECK(hipMemcpyToSymbol(HIP_SYMBOL(DC),C,sizeof C)); HIP_CHECK(hipMemcpyToSymbol(HIP_SYMBOL(DB),B,sizeof B));
  unsigned long long *out=nullptr,zero=0,result=0; HIP_CHECK(hipMalloc(&out,sizeof *out)); HIP_CHECK(hipMemcpy(out,&zero,sizeof zero,hipMemcpyHostToDevice));
  hipLaunchKernelGGL(score,dim3(4096),dim3(256),0,0,out); hipError_t err=hipDeviceSynchronize();
  if(err!=hipSuccess){std::cerr<<hipGetErrorString(err)<<"\n";return 1;} HIP_CHECK(hipMemcpy(&result,out,sizeof result,hipMemcpyDeviceToHost)); HIP_CHECK(hipFree(out));
  unsigned x=(unsigned)(result&(XS-1)),y=(unsigned)((result>>13)&(YS-1)); int energy=host_energy(C,B,x,y);
  std::cout<<"{\"n\":28,\"source_phi\":21,\"phi\":"<<(result>>32)<<",\"energy\":"<<energy<<",\"xmask\":"<<x<<",\"ymask\":"<<y<<",\"workers\":\"hip\"}\n";
}
