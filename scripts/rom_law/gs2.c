// Exact max/min of x^T J x over {+-1}^n with x_{n-1}=+1; optional split: x_{n-2}=+1 (t=0) or -1 (t=1).
// usage: gs2 file [t]   (no t: full enumeration of 2^(n-1) states)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"rb"); int n; if(fread(&n,4,1,f)!=1) return 1;
  double*J=malloc(sizeof(double)*n*n); if(fread(J,8,n*n,f)!=(size_t)(n*n)) return 1; fclose(f);
  int split=(argc>2), t=split?atoi(argv[2]):0;
  double*h=malloc(sizeof(double)*n); double*x=malloc(sizeof(double)*n);
  double tr=0; for(int i=0;i<n;i++){x[i]=1.0; tr+=J[i*n+i];}
  if(split&&t) x[n-2]=-1.0;
  for(int i=0;i<n;i++){double s=0; for(int j=0;j<n;j++) if(j!=i) s+=J[i*n+j]*x[j]; h[i]=s;}
  double E=0; for(int i=0;i<n;i++) E+=x[i]*h[i];
  double mx=E,mn=E; uint64_t total=1ULL<<(split?n-2:n-1);
  for(uint64_t c=1;c<total;c++){
    int k=__builtin_ctzll(c); double xk=x[k];
    E-=4.0*xk*h[k];
    const double*Jk=J+(size_t)k*n; double tt=2.0*xk;
    for(int j=0;j<n;j++) h[j]-=tt*Jk[j];
    h[k]+=tt*Jk[k]; x[k]=-xk;
    if(E>mx)mx=E; if(E<mn)mn=E;
  }
  printf("%.6f %.6f\n",mx+tr,mn+tr); return 0;}
