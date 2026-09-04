#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <vector>
static int inv11(int a){a=(a%11+11)%11;for(int b=1;b<11;b++)if(a*b%11==1)return b;return 0;}
static int direction(int u,int v){int x=u/11,y=u%11,X=v/11,Y=v%11;int dx=(X-x+11)%11,dy=(Y-y+11)%11;if(!dx)return 11;return dy*inv11(dx)%11;}
static std::array<int,11> line(int d,int c){std::array<int,11>L{};if(d==11)for(int y=0;y<11;y++)L[y]=11*c+y;else for(int x=0;x<11;x++)L[x]=11*x+(d*x+c)%11;return L;}
struct Set4{std::array<unsigned char,4>x;std::array<uint64_t,4>sig;bool has0;};
static uint32_t block(const Set4&s,int b){int lo=121*b/10,hi=121*(b+1)/10;uint32_t z=0;for(int v=lo;v<hi;v++){int bit=2*v;int val=(s.sig[bit>>6]>>(bit&63))&3ULL;z|=val<<(2*(v-lo));}return z;}
int main(int argc,char**argv){
 int only=argc>1?std::stoi(argv[1]):-1;
 const int reps[4]={63,95,111,119};
 std::array<std::array<int,11>,132>L{};for(int d=0;d<12;d++)for(int c=0;c<11;c++)L[11*d+c]=line(d,c);
 for(int mask:reps){
  if(only>=0&&mask!=only)continue;
  std::array<std::array<unsigned char,121>,121>A{};for(int u=0;u<121;u++)for(int v=u+1;v<121;v++)A[u][v]=A[v][u]=(mask>>direction(u,v))&1;
  std::vector<Set4>S;S.reserve(43560);
  for(int li=0;li<132;li++)for(int a=0;a<8;a++)for(int b=a+1;b<9;b++)for(int c=b+1;c<10;c++)for(int d=c+1;d<11;d++){
   Set4 s{{(unsigned char)L[li][a],(unsigned char)L[li][b],(unsigned char)L[li][c],(unsigned char)L[li][d]}, {}, false};
   s.has0=std::find(s.x.begin(),s.x.end(),0)!=s.x.end();
   for(int v=0;v<121;v++){int q=0;for(int x:s.x)q+=A[v][x];if(q==4)q=0;int bit=2*v;s.sig[bit>>6]|=(uint64_t)q<<(bit&63);}
   S.push_back(s);
  }
  std::array<std::vector<std::pair<uint64_t,int>>,45> idx;int pp=0;
  for(int bi=0;bi<10;bi++)for(int bj=bi+1;bj<10;bj++){
   auto&V=idx[pp++];V.reserve(S.size());for(int z=0;z<(int)S.size();z++){uint64_t key=((uint64_t)block(S[z],bi)<<32)|block(S[z],bj);V.push_back({key,z});}std::sort(V.begin(),V.end());
  }
  std::unordered_set<uint64_t> seen;long long exact=0;
  for(int i=0;i<(int)S.size();i++)if(S[i].has0){
   pp=0;for(int bi=0;bi<10;bi++)for(int bj=bi+1;bj<10;bj++){
    uint64_t key=((uint64_t)block(S[i],bi)<<32)|block(S[i],bj);auto&V=idx[pp++];
    auto lo=std::lower_bound(V.begin(),V.end(),std::pair<uint64_t,int>{key,-1});auto hi=std::upper_bound(V.begin(),V.end(),std::pair<uint64_t,int>{key,INT32_MAX});
    for(auto it=lo;it!=hi;it++){int j=it->second;if(j==i)continue;uint64_t pair=((uint64_t)std::min(i,j)<<32)|std::max(i,j);if(!seen.insert(pair).second)continue;exact++;
     std::array<char,121> in{};bool ok=true;for(int x:S[i].x)in[x]=1;for(int x:S[j].x){if(in[x])ok=false;in[x]=2;}if(!ok)continue;
     int deg1=0,deg2=0;for(int y:S[i].x)deg1+=A[S[i].x[0]][y];for(int y:S[j].x)deg2+=A[S[j].x[0]][y];if(deg1!=deg2)continue;
     for(int x:S[i].x){int q=0;for(int y:S[j].x)q+=A[x][y];if(x==S[i].x[0])deg1=q;else if(q!=deg1)ok=false;}for(int x:S[j].x){int q=0;for(int y:S[i].x)q+=A[x][y];if(q!=deg1)ok=false;}if(!ok)continue;
     std::vector<int>D1,D2;for(int v=0;v<121&&ok;v++)if(!in[v]){int a=0,b=0;for(int x:S[i].x)a+=A[v][x];for(int x:S[j].x)b+=A[v][x];if(a==b)continue;if(a==4&&b==0)D1.push_back(v);else if(a==0&&b==4)D2.push_back(v);else ok=false;}if(!ok||D1.empty()&&D2.empty())continue;
     auto B=A;for(int v:D1){for(int x:S[i].x)B[v][x]=B[x][v]=0;for(int x:S[j].x)B[v][x]=B[x][v]=1;}for(int v:D2){for(int x:S[j].x)B[v][x]=B[x][v]=0;for(int x:S[i].x)B[v][x]=B[x][v]=1;}
     for(int u=0;u<121&&ok;u++){int d=0;for(int v=0;v<121;v++)d+=B[u][v];if(d!=60)ok=false;}for(int u=0;u<121&&ok;u++)for(int v=0;v<121&&ok;v++){int q=0;for(int k=0;k<121;k++)q+=B[u][k]*B[k][v];int w=u==v?60:(B[u][v]?29:30);if(q!=w)ok=false;}if(!ok)continue;
     std::cout<<"FOUND mask="<<mask<<" C1";for(int x:S[i].x)std::cout<<" "<<x;std::cout<<" C2";for(int x:S[j].x)std::cout<<" "<<x;std::cout<<" D="<<D1.size()<<","<<D2.size()<<"\n";
     std::ofstream f("/tmp/pn_wqh4_lines_candidate.csv");for(int u=0;u<121;u++){for(int v=0;v<121;v++){if(v)f<<',';f<<(int)B[u][v];}f<<'\n';}return 0;
    }
   }
  }
  std::cerr<<"mask="<<mask<<" exact_candidates="<<exact<<"\n";
 }
 std::cout<<"NONE\n";
}
