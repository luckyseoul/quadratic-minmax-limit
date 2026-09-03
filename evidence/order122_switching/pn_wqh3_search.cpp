#include <algorithm>
#include <array>
#include <bit>
#include <fstream>
#include <iostream>
#include <vector>
static int inv11(int a){a=(a%11+11)%11;for(int b=1;b<11;b++)if(a*b%11==1)return b;return 0;}
static int direction(int u,int v){int x=u/11,y=u%11,X=v/11,Y=v%11;int dx=(X-x+11)%11,dy=(Y-y+11)%11;if(!dx)return 11;return dy*inv11(dx)%11;}
int main(int argc,char**argv){
 int only=argc>1?std::stoi(argv[1]):-1;long long total=0;
 for(int mask=0;mask<4096;mask++)if(std::popcount((unsigned)mask)==6&&(only<0||mask==only)){
  std::array<std::array<unsigned char,121>,121>A{};for(int u=0;u<121;u++)for(int v=u+1;v<121;v++)A[u][v]=A[v][u]=(mask>>direction(u,v))&1;
  long long cand=0;int goodc1=0;
  for(int b=1;b<120;b++)for(int c=b+1;c<121;c++){
   int C1[3]={0,b,c};int e1=A[0][b]+A[0][c]+A[b][c];if(e1!=0&&e1!=3)continue;goodc1++;
   std::array<std::vector<int>,8> cls;
   for(int v=1;v<121;v++)if(v!=b&&v!=c){int p=A[v][0]|(A[v][b]<<1)|(A[v][c]<<2);cls[p].push_back(v);}
   auto tryc=[&](int x,int y,int z)->bool{
    if(x==y||x==z||y==z)return false; if(x>y)std::swap(x,y);if(y>z)std::swap(y,z);if(x>y)std::swap(x,y);
    int e2=A[x][y]+A[x][z]+A[y][z];if(e2!=e1)return false;
    cand++;int D1=0,D2=0;bool ok=true;
    for(int v=1;v<121&&ok;v++)if(v!=b&&v!=c&&v!=x&&v!=y&&v!=z){int a=A[v][0]+A[v][b]+A[v][c],d=A[v][x]+A[v][y]+A[v][z];if(a==d)continue;if(a==3&&d==0){D1++;continue;}if(a==0&&d==3){D2++;continue;}ok=false;}
    if(!ok||D1+D2==0)return false;
    auto B=A;int C2[3]={x,y,z};
    for(int v=1;v<121;v++)if(v!=b&&v!=c&&v!=x&&v!=y&&v!=z){int a=A[v][0]+A[v][b]+A[v][c],d=A[v][x]+A[v][y]+A[v][z];if((a==3&&d==0)||(a==0&&d==3))for(int i=0;i<3;i++){B[v][C1[i]]=B[C1[i]][v]=!B[v][C1[i]];B[v][C2[i]]=B[C2[i]][v]=!B[v][C2[i]];}}
    for(int u=0;u<121&&ok;u++){int deg=0;for(int v=0;v<121;v++)deg+=B[u][v];if(deg!=60)ok=false;}
    for(int u=0;u<121&&ok;u++)for(int v=0;v<121&&ok;v++){int q=0;for(int k=0;k<121;k++)q+=B[u][k]*B[k][v];int w=u==v?60:(B[u][v]?29:30);if(q!=w)ok=false;}
    if(!ok){std::cerr<<"BUG cospectral but not SRG\n";return false;}
    std::cout<<"FOUND mask="<<mask<<" C1=0,"<<b<<","<<c<<" C2="<<x<<","<<y<<","<<z<<" D="<<D1<<","<<D2<<"\n";
    std::ofstream f("/tmp/pn_wqh3_candidate.csv");for(int u=0;u<121;u++){for(int v=0;v<121;v++){if(v)f<<',';f<<(int)B[u][v];}f<<'\n';}return true;
   };
   // bipartite degree 0 or 3
   for(int p: {0,7}){auto &V=cls[p];for(int i=0;i<(int)V.size();i++)for(int j=i+1;j<(int)V.size();j++)for(int k=j+1;k<(int)V.size();k++)if(tryc(V[i],V[j],V[k]))return 0;}
   // bipartite degree 1 or 2: one vertex from each required pattern.
   for(auto ps: {std::array<int,3>{1,2,4},std::array<int,3>{6,5,3}})
    for(int x:cls[ps[0]])for(int y:cls[ps[1]])for(int z:cls[ps[2]])if(tryc(x,y,z))return 0;
  }
  total+=cand;std::cerr<<"mask="<<mask<<" C1="<<goodc1<<" cand="<<cand<<" total="<<total<<"\n";
 }
 std::cout<<"NONE candidates="<<total<<"\n";
}

