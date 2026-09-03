#include <algorithm>
#include <array>
#include <bit>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

static int inv11(int a){a=(a%11+11)%11;for(int b=1;b<11;b++)if(a*b%11==1)return b;return 0;}
static int direction(int u,int v){int x=u/11,y=u%11,X=v/11,Y=v%11;int dx=(X-x+11)%11,dy=(Y-y+11)%11;if(!dx)return 11;return dy*inv11(dx)%11;}
static std::array<int,11> line(int d,int c){
  std::array<int,11> L{};
  if(d==11){for(int y=0;y<11;y++)L[y]=11*c+y;}
  else {for(int x=0;x<11;x++)L[x]=11*x+((d*x+c)%11);}
  return L;
}
int main(){
  std::array<std::array<int,11>,132> L{};
  for(int d=0;d<12;d++)for(int c=0;c<11;c++)L[11*d+c]=line(d,c);
  long long tested=0,wqh=0;
  for(int mask=0;mask<(1<<12);mask++){
    if(std::popcount((unsigned)mask)!=6)continue;
    std::array<std::array<unsigned char,121>,121>A{};
    for(int u=0;u<121;u++)for(int v=u+1;v<121;v++)A[u][v]=A[v][u]=(mask>>direction(u,v))&1;
    for(int i=0;i<132;i++)for(int j=i+1;j<132;j++){
      tested++;
      std::vector<int>C1(L[i].begin(),L[i].end()),C2(L[j].begin(),L[j].end());
      std::set<int>s1(C1.begin(),C1.end()),s2(C2.begin(),C2.end());
      std::vector<int> inter; std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::back_inserter(inter));
      if(inter.size()==1){int p=inter[0];C1.erase(std::find(C1.begin(),C1.end(),p));C2.erase(std::find(C2.begin(),C2.end(),p));}
      else if(inter.size()!=0)continue;
      int q=C1.size(); if((int)C2.size()!=q)continue;
      std::array<char,121> in{};for(int x:C1)in[x]=1;for(int x:C2)in[x]=2;
      auto cnt=[&](int u,const std::vector<int>&C){int z=0;for(int v:C)z+=A[u][v];return z;};
      int target=cnt(C1[0],C1)-cnt(C1[0],C2); bool ok=true;
      for(int u:C1)if(cnt(u,C1)-cnt(u,C2)!=target)ok=false;
      for(int v:C2)if(cnt(v,C2)-cnt(v,C1)!=target)ok=false;
      if(!ok)continue;
      std::vector<int>D1,D2;
      for(int v=0;v<121&&ok;v++)if(!in[v]){
        int a=cnt(v,C1),b=cnt(v,C2);
        if(a==q&&b==0)D1.push_back(v);
        else if(a==0&&b==q)D2.push_back(v);
        else if(a!=b)ok=false;
      }
      if(!ok||D1.empty()||D2.empty())continue;
      wqh++;
      auto B=A;
      for(int v:D1){for(int x:C1)B[v][x]=B[x][v]=0;for(int x:C2)B[v][x]=B[x][v]=1;}
      for(int v:D2){for(int x:C2)B[v][x]=B[x][v]=0;for(int x:C1)B[v][x]=B[x][v]=1;}
      for(int u=0;u<121&&ok;u++){int deg=0;for(int v=0;v<121;v++)deg+=B[u][v];if(deg!=60)ok=false;}
      if(!ok)continue;
      // Exact conference-SRG identity.
      for(int u=0;u<121&&ok;u++)for(int v=0;v<121&&ok;v++){
        int z=0;for(int k=0;k<121;k++)z+=B[u][k]*B[k][v];
        int want=(u==v)?60:(B[u][v]?29:30);if(z!=want)ok=false;
      }
      if(!ok)continue;
      std::cout<<"FOUND mask="<<mask<<" lines="<<i<<","<<j<<" q="<<q<<" D="<<D1.size()<<","<<D2.size()<<"\n";
      std::cout<<"C1";for(int x:C1)std::cout<<" "<<x;std::cout<<"\nC2";for(int x:C2)std::cout<<" "<<x;std::cout<<"\nD1";for(int x:D1)std::cout<<" "<<x;std::cout<<"\nD2";for(int x:D2)std::cout<<" "<<x;std::cout<<"\n";
      std::ofstream f("/tmp/pn_wqh_candidate.csv");
      for(int u=0;u<121;u++){for(int v=0;v<121;v++){if(v)f<<',';f<<(int)B[u][v];}f<<'\n';}
      return 0;
    }
  }
  std::cout<<"NONE tested="<<tested<<" wqh_nontrivial="<<wqh<<"\n";
}

