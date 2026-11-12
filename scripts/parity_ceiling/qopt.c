// qopt: optimise V(e) = sum_{a<b} R_ab e_a e_b over e in {+-1}^m for an integer symmetric matrix R (diagonal ignored).
// usage: qopt exact file [part nparts]   -> exact max and min (Gray code, e_{m-1}=+1 fixed)
//        qopt sa file restarts sweeps seed -> annealing lower bounds for max and for -min
// file: m then m*m integers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
static uint64_t s[2];
static inline uint64_t rnd(void){uint64_t s1=s[0],s0=s[1];s[0]=s0;s1^=s1<<23;s[1]=s1^s0^(s1>>17)^(s0>>26);return s[1]+s0;}
static inline double urand(void){return (rnd()>>11)*(1.0/9007199254740992.0);}
int main(int argc,char**argv){
  FILE*f=fopen(argv[2],"r"); int m; if(fscanf(f,"%d",&m)!=1) return 1;
  long *R=malloc(sizeof(long)*m*m); for(int i=0;i<m*m;i++) if(fscanf(f,"%ld",&R[i])!=1) return 1; fclose(f);
  for(int i=0;i<m;i++) R[i*m+i]=0;
  if(!strcmp(argv[1],"exact")){
    int part=argc>3?atoi(argv[3]):0, nparts=argc>4?atoi(argv[4]):1; // nparts power of two, fixes top free bits
    int kb=0; while((1<<kb)<nparts) kb++;
    int free=m-1-kb; // Gray bits
    long *e=malloc(sizeof(long)*m), *F=malloc(sizeof(long)*m);
    for(int i=0;i<m;i++) e[i]=1;
    for(int b=0;b<kb;b++) if((part>>b)&1) e[m-2-b]=-1;
    long V=0; for(int a=0;a<m;a++){long fa=0; for(int b=0;b<m;b++) fa+=R[a*m+b]*e[b]; F[a]=fa; V+=e[a]*fa;} V/=2;
    long best=V,worst=V; uint64_t bestg=0,worstg=0;
    uint64_t N=1ULL<<free;
    for(uint64_t g=1;g<N;g++){
      int k=__builtin_ctzll(g);
      long d=e[k]*F[k]; V-=2*d; e[k]=-e[k]; long t=2*e[k]; const long*row=R+(size_t)k*m;
      for(int b=0;b<m;b++) F[b]+=row[b]*t;
      if(V>best){best=V;bestg=g;} if(V<worst){worst=V;worstg=g;}
    }
    printf("exact m=%d part=%d/%d max=%ld min=%ld gmax=%llu gmin=%llu\n",m,part,nparts,best,worst,(unsigned long long)bestg,(unsigned long long)worstg);
    return 0;
  }
  int Rn=atoi(argv[3]), SW=atoi(argv[4]); uint64_t seed=strtoull(argv[5],0,10);
  s[0]=seed*0x9E3779B97F4A7C15ULL+1; s[1]=seed^0xD1B54A32D192ED03ULL; for(int i=0;i<20;i++)rnd();
  long *e=malloc(sizeof(long)*m), *F=malloc(sizeof(long)*m); long res[2]; int hits[2];
  double scale=0; for(int i=0;i<m*m;i++) scale+=fabs((double)R[i]); scale/=m; // mean |row sum| scale
  for(int sgn=0;sgn<2;sgn++){
    long best=-(1L<<60); int nb=0; long sg=sgn?-1:1;
    for(int r=0;r<Rn;r++){
      for(int i=0;i<m;i++) e[i]=(rnd()&1)?1:-1;
      long V=0; for(int a=0;a<m;a++){long fa=0; for(int b=0;b<m;b++) fa+=sg*R[a*m+b]*e[b]; F[a]=fa; V+=e[a]*fa;} V/=2;
      double T=0.5*scale, Tmin=0.01*scale, fac=pow(Tmin/T,1.0/(SW-1));
      for(int sw=0;sw<SW;sw++){ for(int q=0;q<m;q++){ int k=(int)(rnd()%m); long d=e[k]*F[k];
          if(d<=0||urand()<exp(-2.0*d/T)){ V-=2*d; e[k]=-e[k]; long t=2*e[k]*sg; const long*row=R+(size_t)k*m; for(int b=0;b<m;b++) F[b]+=row[b]*t; } }
        T*=fac; }
      int imp=1; while(imp){imp=0; for(int k=0;k<m;k++){long d=e[k]*F[k]; if(d<0){V-=2*d; e[k]=-e[k]; long t=2*e[k]*sg; const long*row=R+(size_t)k*m; for(int b=0;b<m;b++) F[b]+=row[b]*t; imp=1;}}}
      if(V>best){best=V;nb=1;} else if(V==best) nb++;
    }
    res[sgn]=best; hits[sgn]=nb;
  }
  printf("sa m=%d max>=%ld (hits %d/%d) min<=%ld (hits %d/%d)\n",m,res[0],hits[0],Rn,-res[1],hits[1],Rn);
  return 0;
}
