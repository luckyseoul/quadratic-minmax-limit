#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Independent exact certificate for one finite fiber only:
//   p=31, b=7, k=11, constant equianharmonic tangent-conic imbalance,
//   with every degree-six and degree-eight channel required to vanish.
//
// It uses no SAT/MIP library.  Alignment leaves total deficit three.  The
// nine maximal AE atoms are disjoint conic triangles, so their multiplicities
// are reconstructed uniquely.  Maximal compact atoms are enumerated by their
// one off-target orbit; the three-defect all-compact case uses an exact 3SUM
// table on the unsupported off-target projection.

using namespace std;
constexpr int P=31, KPAR=11, NT=29, NO=225, NM=7;
using A7=array<int,NM>;
using A20=array<int,20>;
using A29=array<int,NT>;

int mod(int x){x%=P; if(x<0)x+=P; return x;}
int pw(int a,int n){int r=1;a=mod(a);while(n){if(n&1)r=r*a%P;a=a*a%P;n>>=1;}return r;}
pair<int,int> edge(int a,int b){a=mod(a);b=mod(b);assert(a!=b);if(a>b)swap(a,b);return {a,b};}
pair<int,int> negedge(pair<int,int> e){return edge(-e.first,-e.second);}
struct Orb{int id,sign;};
map<pair<int,int>,int> oid;
vector<pair<int,int>> oreps;
Orb orbit(pair<int,int> e){auto ne=negedge(e);if(ne==e)return {-1,0};auto rep=min(e,ne);return {oid.at(rep),e==rep?1:-1};}

struct Atom{
  bool ae;
  array<int,3> tri;
  int distinguished=-1;
  int score=0,deficit=0;
  A29 target{};
  vector<pair<int,int>> off;
  A7 mom{};
};
vector<Atom> aeall, kvals[5]; // K scores -2..2 at index score+2
map<int,int> target_by_oid;
array<int,NO> target_t{};
array<int,NO> target_index{};
array<array<int,3>,9> cycle_ids;
array<A7,9> cycle_mom;
array<const Atom*,9> cycle_atoms{};
A7 mconst{};
array<int,2> broken{};

void add_edge(array<int,NO>& con,int u,int v,int coefficient){
  auto o=orbit(edge(u,v));if(o.id>=0)con[o.id]+=coefficient*o.sign;
}
int edge_moment(int s,int t,int degree,int channel){
  return pw(s-t,2)*pw(s*t,channel)%P*pw(s+t,degree-2-2*channel)%P;
}

Atom make_atom(array<int,3> tri,int distinguished){
  Atom a; a.ae=distinguished<0;a.tri=tri;a.distinguished=distinguished;
  array<int,NO> con{};
  vector<tuple<int,int,int>> raw;
  if(a.ae){
    raw={{tri[0],tri[1],1},{tri[0],tri[2],1},{tri[1],tri[2],1}};
  }else{
    vector<int> q;for(int x:tri)if(x!=distinguished)q.push_back(x);assert(q.size()==2);
    raw={{q[0],q[1],1},{q[0],distinguished,-1},{q[1],distinguished,-1}};
  }
  for(auto [u,v,c]:raw)add_edge(con,u,v,c);
  for(auto [u,v,c]:raw){
    int coord=0;
    for(int d:{6,8})for(int j=0;j<d/2;++j){
      a.mom[coord]=mod(a.mom[coord]+c*edge_moment(u,v,d,j));
      ++coord;
    }
  }
  for(int o=0;o<NO;++o)if(con[o]){
    int ti=target_index[o];
    if(ti>=0){a.target[ti]=con[o]*target_t[o];a.score+=a.target[ti];}
    else a.off.push_back({o,con[o]});
  }
  a.deficit=(a.ae?3:2)-a.score;
  return a;
}

struct Feature{
  A20 lin{}; A7 adj{}; A29 kt{}; array<int,NO> off{};
};
Feature feature(const vector<const Atom*>& rs){
  Feature f;A7 mom{};
  for(auto a:rs){
    for(int i=0;i<NT;++i)f.kt[i]+=a->target[i];
    for(auto [o,x]:a->off)f.off[o]+=x;
    for(int j=0;j<NM;++j)mom[j]=mod(mom[j]+a->mom[j]);
  }
  int z=0;
  for(int c=0;c<9;++c){auto ids=cycle_ids[c];f.lin[z++]=f.kt[ids[1]]-f.kt[ids[0]];f.lin[z++]=f.kt[ids[2]]-f.kt[ids[0]];}
  f.lin[z++]=f.kt[broken[0]];f.lin[z++]=f.kt[broken[1]];assert(z==20);
  f.adj=mom;
  for(int c=0;c<9;++c){int v=f.kt[cycle_ids[c][0]];for(int j=0;j<NM;++j)f.adj[j]=mod(f.adj[j]-v*cycle_mom[c][j]);}
  return f;
}
bool offzero(const Feature&f){for(int x:f.off)if(x)return false;return true;}
bool linok(const Feature&f){for(int i=0;i<18;++i)if(f.lin[i])return false;return f.lin[18]==1&&f.lin[19]==1;}
bool adjok(const Feature&f){for(int j=0;j<NM;++j)if(mod(f.adj[j]+mconst[j]))return false;return true;}

vector<const Atom*> maxnull;
array<array<vector<const Atom*>,2>,NO> maxsigned; // sign - at0, + at1
vector<int> maxorbits;

void choose_multiset(const vector<const Atom*>& pool,int count,int start,vector<const Atom*>& cur,const function<void()>& cb){
  if(count==0){cb();return;}
  for(int i=start;i<(int)pool.size();++i){cur.push_back(pool[i]);choose_multiset(pool,count-1,i,cur,cb);cur.pop_back();}
}

// Exact enumeration of score-2 compact multisets with prescribed off sum.
void max_completions(const array<int,NO>& required,int count,const function<bool(const vector<const Atom*>&)>& cb){
  vector<pair<int,int>> req;
  int forced_count=0;
  for(int o=0;o<NO;++o)if(required[o]){
    int s=required[o]>0?1:-1,n=abs(required[o]);
    if(maxsigned[o][s>0?1:0].empty())return;
    req.push_back({o,required[o]});forced_count+=n;
  }
  if(forced_count>count)return;
  vector<const Atom*> cur;
  function<bool(int)> forced_rec=[&](int qi){
    if(qi<(int)req.size()){
      auto [o,v]=req[qi];int n=abs(v),si=v>0?1:0;
      bool stop=false;
      choose_multiset(maxsigned[o][si],n,0,cur,[&](){if(!stop)stop=forced_rec(qi+1);});
      return stop;
    }
    int rem=count-forced_count;
    for(int np=0;np<=rem/2;++np){
      int nn=rem-2*np;
      vector<int> chosen_orbits;
      function<bool(int,int)> pair_orbit_rec=[&](int left,int start){
        if(left){
          for(int i=start;i<(int)maxorbits.size();++i){chosen_orbits.push_back(maxorbits[i]);if(pair_orbit_rec(left-1,i))return true;chosen_orbits.pop_back();}
          return false;
        }
        map<int,int> counts;for(int o:chosen_orbits)counts[o]++;
        vector<pair<int,int>> groups(counts.begin(),counts.end());
        function<bool(int)> pair_choice_rec=[&](int gi){
          if(gi<(int)groups.size()){
            auto [o,n]=groups[gi];bool stop=false;
            choose_multiset(maxsigned[o][0],n,0,cur,[&](){
              choose_multiset(maxsigned[o][1],n,0,cur,[&](){if(!stop)stop=pair_choice_rec(gi+1);});
            });return stop;
          }
          bool stop=false;choose_multiset(maxnull,nn,0,cur,[&](){if(!stop)stop=cb(cur);});return stop;
        };
        return pair_choice_rec(0);
      };
      if(pair_orbit_rec(np,0))return true;
    }
    return false;
  };
  forced_rec(0);
}

void print_atom(const Atom* a){
  cout<<(a->ae?"AE":"K")<<" ("<<a->tri[0]<<","<<a->tri[1]<<","<<a->tri[2]<<")";
  if(!a->ae)cout<<" d="<<a->distinguished;
  cout<<" score="<<a->score<<"\n";
}

bool check_hit(const vector<const Atom*>& exceptional,const vector<const Atom*>& maximal,uint64_t tested,const string&part){
  vector<const Atom*> selected=exceptional;selected.insert(selected.end(),maximal.begin(),maximal.end());
  auto f=feature(selected);if(!offzero(f)||!linok(f)||!adjok(f))return false;
  int ea=0;for(auto a:exceptional)ea+=a->ae;
  int asum=0;
  for(int c=0;c<9;++c){int m=1-f.kt[cycle_ids[c][0]];if(m<0)return false;asum+=m;}
  if(asum!=6-ea)return false;
  cout<<"SAT partition="<<part<<" tested="<<tested<<"\nEXCEPTIONAL\n";
  for(auto a:exceptional)print_atom(a);
  cout<<"MAXIMAL K\n";
  for(auto a:maximal)print_atom(a);
  cout<<"AE CYCLE MULTIPLICITIES\n";
  for(int c=0;c<9;++c){int m=1-f.kt[cycle_ids[c][0]];print_atom(cycle_atoms[c]);cout<<"multiplicity="<<m<<"\n";}
  return true;
}

void init(){
  target_index.fill(-1);
  set<pair<int,int>> reps;
  for(int a=0;a<P;++a)for(int b=a+1;b<P;++b){auto e=edge(a,b),ne=negedge(e);if(e!=ne)reps.insert(min(e,ne));}
  oreps.assign(reps.begin(),reps.end());assert(oreps.size()==NO);for(int i=0;i<NO;++i)oid[oreps[i]]=i;
  int inv2=pw(2,P-2);
  for(int z=0;z<P;++z)if(z!=0&&z!=1){int u=mod(((1+KPAR)*z-KPAR)*inv2),v=mod(((1-KPAR)*z+KPAR)*inv2);auto o=orbit(edge(u,v));target_by_oid[o.id]=o.sign;target_t[o.id]=o.sign;}
  assert(target_by_oid.size()==NT);int ti=0;for(auto [o,s]:target_by_oid)target_index[o]=ti++;
  for(int a=0;a<P;++a)for(int b=a+1;b<P;++b)for(int c=b+1;c<P;++c){
    array<int,3> t{a,b,c};aeall.push_back(make_atom(t,-1));
    for(int d:t){auto ka=make_atom(t,d);assert(-2<=ka.score&&ka.score<=2);kvals[ka.score+2].push_back(std::move(ka));}
  }
  vector<const Atom*> cycles;
  array<int,7> ascores{};for(auto &a:aeall){assert(-3<=a.score&&a.score<=3);ascores[a.score+3]++;if(a.score==3)cycles.push_back(&a);}
  assert((ascores==array<int,7>{9,1,702,3071,702,1,9}));
  array<int,5> kscores{};for(int i=0;i<5;++i)kscores[i]=kvals[i].size();
  assert((kscores==array<int,5>{111,2133,8997,2133,111}));
  assert(cycles.size()==9);
  set<int> used;
  for(int c=0;c<9;++c){cycle_atoms[c]=cycles[c];int q=0;for(int i=0;i<NT;++i)if(cycles[c]->target[i]){cycle_ids[c][q++]=i;used.insert(i);}assert(q==3);cycle_mom[c]=cycles[c]->mom;for(int j=0;j<NM;++j)mconst[j]=mod(mconst[j]+cycle_mom[c][j]);}
  int q=0;for(int i=0;i<NT;++i)if(!used.count(i))broken[q++]=i;assert(q==2);
  for(auto &a:kvals[4]){if(a.off.empty())maxnull.push_back(&a);else{assert(a.off.size()==1&&abs(a.off[0].second)==1);int o=a.off[0].first,s=a.off[0].second;maxsigned[o][s>0?1:0].push_back(&a);}}
  for(int o=0;o<NO;++o)if(!maxsigned[o][0].empty()||!maxsigned[o][1].empty()){assert(!maxsigned[o][0].empty()&&!maxsigned[o][1].empty());maxorbits.push_back(o);}
  cerr<<"init cycles="<<cycles.size()<<" maxK="<<kvals[4].size()<<" null="<<maxnull.size()<<" maxorbits="<<maxorbits.size()<<"\n";
}

int main(int argc,char**argv){
  init();string mode=argc>1?argv[1]:"d3";
  if(mode=="manifest"){
    cout<<"schema=p31_equianharmonic_zero68_mitm_manifest_v1\n"
        <<"scope=p31_b7_k11_constant_conic_fiber_only\n"
        <<"status=UNSAT\n"
        <<"all_completions_tested=230314710\n"
        <<"all_edge_hits=17076\n"
        <<"all_zero68_hits=0\n"
        <<"d3_completions=13528344 d3_edge_hits=60\n"
        <<"d2d1_pairs=20697666 d2d1_compatible=79918 d2d1_completions=87840508 d2d1_edge_hits=2160\n"
        <<"d1_with_ae_exceptions=2278045 d1_with_ae_compatible=24828 d1_with_ae_completions=20465801 d1_with_ae_edge_hits=392\n"
        <<"d1_all_k_multisets=1619689995 d1_all_k_3sum=2027542 d1_all_k_compatible=1089526 d1_all_k_completions=108480057 d1_all_k_edge_hits=14464\n";
  }else if(mode=="d3"){
    vector<const Atom*> cand;for(auto&a:aeall)if(a.deficit==3)cand.push_back(&a);for(auto&a:kvals[1])cand.push_back(&a); // score -1 => index1
    uint64_t tested=0,linhits=0;
    for(auto ex:cand){
      array<int,NO> req{};for(auto [o,x]:ex->off)req[o]-=x;int nk=7-(!ex->ae);
      vector<const Atom*> exceptional{ex};
      bool hit=false;
      max_completions(req,nk,[&](const vector<const Atom*>&mx){++tested;auto f=feature([&](){vector<const Atom*>s=exceptional;s.insert(s.end(),mx.begin(),mx.end());return s;}());if(offzero(f)&&linok(f))++linhits;if(check_hit(exceptional,mx,tested,"one-d3")){hit=true;return true;}return false;});
      if(hit)return 0;
    }
    cout<<"UNSAT_IN_PARTITION one-d3 tested="<<tested<<" linhits="<<linhits<<" candidates="<<cand.size()<<"\n";
  }else if(mode=="d2d1"){
    vector<const Atom*> d2,d1;
    for(auto&a:aeall)if(a.deficit==2)d2.push_back(&a);
    for(auto&a:kvals[2])d2.push_back(&a); // K score 0
    for(auto&a:aeall)if(a.deficit==1)d1.push_back(&a);
    for(auto&a:kvals[3])d1.push_back(&a); // K score 1
    int lo=argc>2?stoi(argv[2]):0,hi=argc>3?stoi(argv[3]):(int)d2.size();
    lo=max(0,lo);hi=min((int)d2.size(),hi);
    uint64_t pairs=0,compatible=0,tested=0,linhits=0;
    array<int,NO> req{};vector<int> touched;
    for(int ii=lo;ii<hi;++ii){auto x=d2[ii];
      for(auto y:d1){++pairs;touched.clear();
        auto addreq=[&](const Atom*a){for(auto [o,v]:a->off){if(req[o]==0)touched.push_back(o);req[o]-=v;}};
        addreq(x);addreq(y);int nk=7-(!x->ae)-(!y->ae),forced=0;bool ok=true;
        for(int o:touched)if(req[o]){forced+=abs(req[o]);if(maxsigned[o][req[o]>0?1:0].empty())ok=false;}
        if(forced>nk)ok=false;
        if(ok){++compatible;vector<const Atom*> exceptional{x,y};bool hit=false;
          max_completions(req,nk,[&](const vector<const Atom*>&mx){++tested;vector<const Atom*> all=exceptional;all.insert(all.end(),mx.begin(),mx.end());auto f=feature(all);if(offzero(f)&&linok(f))++linhits;if(check_hit(exceptional,mx,tested,"d2+d1")){hit=true;return true;}return false;});
          if(hit)return 0;
        }
        for(int o:touched)req[o]=0;
      }
      if((ii-lo)%250==0)cerr<<"progress d2="<<ii<<"/"<<hi<<" pairs="<<pairs<<" compatible="<<compatible<<" tested="<<tested<<"\n";
    }
    cout<<"UNSAT_IN_PARTITION d2+d1 range="<<lo<<":"<<hi<<" d2total="<<d2.size()<<" d1total="<<d1.size()<<" pairs="<<pairs<<" compatible="<<compatible<<" tested="<<tested<<" linhits="<<linhits<<"\n";
  }else if(mode=="d1ae"){
    vector<const Atom*> ak;
    for(auto&a:aeall)if(a.deficit==1)ak.push_back(&a);
    assert(ak.size()==1);const Atom* A=ak[0];auto &ks=kvals[3];
    int lo=argc>2?stoi(argv[2]):0,hi=argc>3?stoi(argv[3]):(int)ks.size();lo=max(0,lo);hi=min((int)ks.size(),hi);
    uint64_t exceptions=0,compatible=0,tested=0,linhits=0;
    auto run_exception=[&](vector<const Atom*> ex){
      ++exceptions;array<int,NO> req{};int nk=7;for(auto a:ex){nk-=!a->ae;for(auto[o,v]:a->off)req[o]-=v;}
      int forced=0;bool ok=true;for(int o=0;o<NO;++o)if(req[o]){forced+=abs(req[o]);if(maxsigned[o][req[o]>0?1:0].empty())ok=false;}
      if(!ok||forced>nk)return false;
      ++compatible;bool hit=false;
      max_completions(req,nk,[&](const vector<const Atom*>&mx){++tested;vector<const Atom*> all=ex;all.insert(all.end(),mx.begin(),mx.end());auto f=feature(all);if(offzero(f)&&linok(f))++linhits;if(check_hit(ex,mx,tested,"three-d1-with-AE")){hit=true;return true;}return false;});return hit;
    };
    // The all-AE and two-AE cases are assigned only to the shard beginning 0.
    if(lo==0){if(run_exception({A,A,A}))return 0;for(auto &k:ks)if(run_exception({A,A,&k}))return 0;}
    for(int i=lo;i<hi;++i){for(int j=i;j<(int)ks.size();++j)if(run_exception({A,&ks[i],&ks[j]}))return 0;if((i-lo)%50==0)cerr<<"progress d1ae i="<<i<<"/"<<hi<<" exceptions="<<exceptions<<" tested="<<tested<<"\n";}
    cout<<"UNSAT_IN_PARTITION three-d1-with-AE range="<<lo<<":"<<hi<<" ktotal="<<ks.size()<<" exceptions="<<exceptions<<" compatible="<<compatible<<" tested="<<tested<<" linhits="<<linhits<<"\n";
  }else if(mode=="d1k"){
    auto &ks=kvals[3];const int n=ks.size();assert(n==2133);
    array<char,NO> supported{};for(int o:maxorbits)supported[o]=1;
    vector<vector<pair<uint8_t,int8_t>>> uoff(n);
    for(int i=0;i<n;++i)for(auto[o,v]:ks[i].off)if(!supported[o])uoff[i].push_back({(uint8_t)o,(int8_t)v});
    auto pairkey=[&](int i,int j){
      string s;auto&a=uoff[i];auto&b=uoff[j];size_t x=0,y=0;
      while(x<a.size()||y<b.size()){
        int o,v;if(y==b.size()||(x<a.size()&&a[x].first<b[y].first)){o=a[x].first;v=a[x++].second;}
        else if(x==a.size()||b[y].first<a[x].first){o=b[y].first;v=b[y++].second;}
        else{o=a[x].first;v=a[x++].second+b[y++].second;}
        if(v){s.push_back((char)o);s.push_back((char)(v+4));}
      }return s;
    };
    auto negkey=[&](int i){string s;for(auto[o,v]:uoff[i]){s.push_back((char)o);s.push_back((char)(-v+4));}return s;};
    unordered_map<string,vector<uint32_t>> pairs;pairs.reserve(1800000);
    for(int i=0;i<n;++i)for(int j=i;j<n;++j)pairs[pairkey(i,j)].push_back(((uint32_t)i<<12)|(uint32_t)j);
    cerr<<"d1k pairtable keys="<<pairs.size()<<" pairs="<<((uint64_t)n*(n+1)/2)<<"\n";
    int lo=argc>2?stoi(argv[2]):0,hi=argc>3?stoi(argv[3]):n;lo=max(0,lo);hi=min(n,hi);
    uint64_t triples=0,compatible=0,tested=0,linhits=0;
    array<int,NO> req{};vector<int> touched;
    for(int k=lo;k<hi;++k){auto it=pairs.find(negkey(k));if(it==pairs.end())continue;
      for(uint32_t packed:it->second){int i=packed>>12,j=packed&4095;if(j>k)continue;++triples;touched.clear();
        for(const Atom*a:{&ks[i],&ks[j],&ks[k]})for(auto[o,v]:a->off){if(req[o]==0)touched.push_back(o);req[o]-=v;}
        sort(touched.begin(),touched.end());touched.erase(unique(touched.begin(),touched.end()),touched.end());
        int forced=0;bool ok=true;for(int o:touched)if(req[o]){forced+=abs(req[o]);if(maxsigned[o][req[o]>0?1:0].empty())ok=false;}
        if(forced>4)ok=false;
        if(ok){++compatible;vector<const Atom*> ex{&ks[i],&ks[j],&ks[k]};bool hit=false;
          max_completions(req,4,[&](const vector<const Atom*>&mx){++tested;vector<const Atom*> all=ex;all.insert(all.end(),mx.begin(),mx.end());auto f=feature(all);if(offzero(f)&&linok(f))++linhits;if(check_hit(ex,mx,tested,"three-K-d1")){hit=true;return true;}return false;});if(hit)return 0;
        }
        for(int o:touched)req[o]=0;
      }
      if((k-lo)%50==0)cerr<<"progress d1k k="<<k<<"/"<<hi<<" triples="<<triples<<" compatible="<<compatible<<" tested="<<tested<<"\n";
    }
    cout<<"UNSAT_IN_PARTITION three-K-d1 range="<<lo<<":"<<hi<<" ktotal="<<n<<" pairkeys="<<pairs.size()<<" triples_after_unsupported_3sum="<<triples<<" compatible="<<compatible<<" tested="<<tested<<" linhits="<<linhits<<"\n";
  }else{
    cerr<<"usage: "<<argv[0]<<" {manifest|d3|d2d1 lo hi|d1ae lo hi|d1k lo hi}\n";
    return 2;
  }
}
