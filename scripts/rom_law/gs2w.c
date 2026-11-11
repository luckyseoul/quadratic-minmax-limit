// Shard walker printing the witness bitmask of the max and min.
// usage: gs2w file s prefix
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"rb"); int n; if(fread(&n,4,1,f)!=1) return 1;
  double*J=malloc(sizeof(double)*n*n); if(fread(J,8,n*n,f)!=(size_t)(n*n)) return 1; fclose(f);
  int s=atoi(argv[2]); uint64_t prefix=(uint64_t)strtoull(argv[3],0,10);
  double*h=malloc(sizeof(double)*n); double*x=malloc(sizeof(double)*n);
  double tr=0; for(int i=0;i<n;i++){x[i]=1.0; tr+=J[i*n+i];}
  for(int i=0;i<s;i++) if((prefix>>i)&1) x[i]=-1.0;
  for(int i=0;i<n;i++){double sm=0; for(int j=0;j<n;j++) if(j!=i) sm+=J[i*n+j]*x[j]; h[i]=sm;}
  double E=0; for(int i=0;i<n;i++) E+=x[i]*h[i];
  double mx=E,mn=E; uint64_t mask=0; for(int i=0;i<n;i++) if(x[i]<0) mask|=1ULL<<i;
  uint64_t mmax=mask, mmin=mask;
  uint64_t total=1ULL<<(n-1-s);
  for(uint64_t c=1;c<total;c++){
    int k=__builtin_ctzll(c)+s; double xk=x[k];
    E-=4.0*xk*h[k];
    const double*Jk=J+(size_t)k*n; double tt=2.0*xk;
    for(int j=0;j<n;j++) h[j]-=tt*Jk[j];
    h[k]+=tt*Jk[k]; x[k]=-xk;
    if(xk<0) mask&=~(1ULL<<k); else mask|=1ULL<<k;
    if(E>mx){mx=E; mmax=mask;}
    if(E<mn){mn=E; mmin=mask;}
  }
  printf("%.6f %llu %.6f %llu\n",mx+tr,(unsigned long long)mmax,mn+tr,(unsigned long long)mmin); return 0;}
