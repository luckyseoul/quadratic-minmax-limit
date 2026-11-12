// Simulated annealing lower bounds for max_x sum_{i<j} C_ij x_i x_j, C = Paley conference matrix of prime order q+1.
// usage: sa q restarts sweeps seed [Tmax Tmin]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
static uint64_t s[2];
static inline uint64_t rnd(void){uint64_t s1=s[0],s0=s[1];s[0]=s0;s1^=s1<<23;s[1]=s1^s0^(s1>>17)^(s0>>26);return s[1]+s0;}
static inline double urand(void){return (rnd()>>11)*(1.0/9007199254740992.0);}
int main(int argc,char**argv){
  int q=atoi(argv[1]); int R=atoi(argv[2]); int SW=atoi(argv[3]); uint64_t seed=strtoull(argv[4],0,10);
  int n=q+1; double Tmax=argc>5?atof(argv[5]):1.2*sqrt((double)q); double Tmin=argc>6?atof(argv[6]):0.35;
  s[0]=seed*0x9E3779B97F4A7C15ULL+1; s[1]=seed^0xD1B54A32D192ED03ULL; for(int i=0;i<20;i++)rnd();
  signed char *C=malloc((size_t)n*n); int *chi=calloc(q,sizeof(int));
  for(int a=1;a<q;a++) chi[(int)(((long long)a*a)%q)]=1;
  for(int a=1;a<q;a++) if(!chi[a]) chi[a]=-1;
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){ int v;
    if(i==j) v=0; else if(i==0||j==0) v=1; else v=chi[((i-j)%q+q)%q]; C[(size_t)i*n+j]=v; }
  int *x=malloc(n*sizeof(int)),*F=malloc(n*sizeof(int)),*bx=malloc(n*sizeof(int)); long best=-1; 
  long *hist=calloc(100000,sizeof(long)); // distribution of final values
  int nb=0;
  for(int r=0;r<R;r++){
    for(int i=0;i<n;i++) x[i]=(rnd()&1)?1:-1;
    long E=0; for(int i=0;i<n;i++){int f=0; const signed char*row=C+(size_t)i*n; for(int j=0;j<n;j++) f+=row[j]*x[j]; F[i]=f; E+=x[i]*f;} E/=2;
    double T=Tmax, fac=pow(Tmin/Tmax,1.0/(SW-1));
    for(int sw=0;sw<SW;sw++){
      for(int k=0;k<n;k++){
        int i=(int)(rnd()%n); int d=x[i]*F[i]; // flipping changes E by -2d
        if(d<=0 || urand()<exp(-2.0*d/T)){
          x[i]=-x[i]; E-=2*d; int xi2=2*x[i]; const signed char*row=C+(size_t)i*n;
          for(int j=0;j<n;j++) F[j]+=row[j]*xi2;
        }
      }
      T*=fac;
    }
    // greedy 1-flip ascent
    int imp=1; while(imp){imp=0; for(int i=0;i<n;i++){int d=x[i]*F[i]; if(d<0){x[i]=-x[i];E-=2*d;int xi2=2*x[i];const signed char*row=C+(size_t)i*n;for(int j=0;j<n;j++)F[j]+=row[j]*xi2;imp=1;}}}
    if(E>=0&&E<100000) hist[E]++;
    if(E>best){best=E; memcpy(bx,x,n*sizeof(int)); nb=1;} else if(E==best) nb++;
  }
  printf("q=%d n=%d best=%ld hits=%d/%d\n",q,n,best,nb,R);
  printf("x="); for(int i=0;i<n;i++) putchar(bx[i]>0?'+':'-'); printf("\n");
  // top of distribution
  printf("dist:"); int shown=0; for(long e=best;e>=0&&shown<6;e--) if(hist[e]){printf(" %ld:%ld",e,hist[e]);shown++;} printf("\n");
  return 0;
}
