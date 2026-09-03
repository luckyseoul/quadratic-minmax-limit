#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using Tri = std::array<int, 3>;
static constexpr int P = 31;
static constexpr int C = 225;

static std::vector<Tri> tris;
static std::vector<std::vector<std::pair<int,int>>> bd;
static std::array<std::array<std::vector<int>,2>, C> cand;
static std::vector<std::array<int,7>> mv;
static std::vector<int> invariant_ids;
static std::vector<int> nonfixed_reps;
static std::array<std::unordered_set<std::uint64_t>,6> fixed_sum_codes;
static std::vector<std::array<int,7>> pair_unit_vectors;
static std::unordered_set<std::uint64_t> pair_unit_codes;
static std::unordered_set<std::uint64_t> pair_pair_codes;
static std::unordered_set<std::uint64_t> pasch_codes;
static std::array<int,C> bal{};
static std::array<int,7> total{};
static std::vector<int> chosen;
static std::array<std::unordered_set<std::string>,7> seen;
static std::uint64_t nodes=0, compact_done=0;
static bool found=false;
static Tri found_compact{};
static std::vector<int> found_ae;
static int found_remaining=0;
static std::array<int,7> found_central_need{};

static Tri norm_tri(int a,int b,int c){
    Tri t{(a%P+P)%P,(b%P+P)%P,(c%P+P)%P};
    std::sort(t.begin(),t.end()); return t;
}
static int pw(int a,int e){ long long r=1,b=(a%P+P)%P; while(e){if(e&1)r=r*b%P;b=b*b%P;e>>=1;}return (int)r; }
static int q(int s,int t,int d,int k){return (long long)pw(s-t,2)*pw(s*t,k)%P*pw(s+t,d-2-2*k)%P;}
static std::uint64_t encode(const std::array<int,7>&v){std::uint64_t z=0,m=1;for(int x:v){z+=m*x;m*=P;}return z;}
static std::array<int,7> decode(std::uint64_t z){std::array<int,7>v{};for(int&i:v){i=z%P;z/=P;}return v;}
static std::array<int,7> addv(const std::array<int,7>&a,const std::array<int,7>&b){std::array<int,7>z{};for(int i=0;i<7;++i)z[i]=(a[i]+b[i])%P;return z;}
static std::array<int,7> subv(const std::array<int,7>&a,const std::array<int,7>&b){std::array<int,7>z{};for(int i=0;i<7;++i)z[i]=(a[i]-b[i]+P)%P;return z;}
static std::array<int,7> ae_moment(const Tri&t){
    std::array<std::pair<int,int>,7> ch{{{6,0},{6,1},{6,2},{8,0},{8,1},{8,2},{8,3}}};
    std::array<int,7> z{}; std::array<std::pair<int,int>,3> es{{{t[0],t[1]},{t[0],t[2]},{t[1],t[2]}}};
    for(int j=0;j<7;++j)for(auto [a,b]:es)z[j]=(z[j]+q(a,b,ch[j].first,ch[j].second))%P;
    return z;
}
static std::array<int,7> compact_moment_vec(const Tri&t){
    std::array<std::pair<int,int>,7> ch{{{6,0},{6,1},{6,2},{8,0},{8,1},{8,2},{8,3}}};
    std::array<int,7> z{};
    for(int j=0;j<7;++j)z[j]=(q(t[0],t[1],ch[j].first,ch[j].second)-q(t[0],t[2],ch[j].first,ch[j].second)-q(t[1],t[2],ch[j].first,ch[j].second)+2*P)%P;
    return z;
}
static std::string key(){
    std::string out;out.reserve(50);
    for(int i=0;i<C;++i)if(bal[i]){out.push_back(char(i));out.push_back(char(bal[i]+16));}
    out.push_back(char(255));
    for(int j=0;j<7;++j)out.push_back(char(total[j]));
    return out;
}

static void finish_invariants(int remaining,int floor_pos=0){
    if(found)return;
    if(remaining==0){
        for(int x:total)if(x%P)return;
        found=true;found_ae=chosen;return;
    }
    for(int pos=floor_pos;pos<(int)invariant_ids.size();++pos){
        int id=invariant_ids[pos];chosen.push_back(id);for(int j=0;j<7;++j)total[j]=(total[j]+mv[id][j])%P;
        finish_invariants(remaining-1,pos);
        for(int j=0;j<7;++j)total[j]=(total[j]-mv[id][j]+P)%P;chosen.pop_back();
        if(found)return;
    }
}

static bool central_completion(int blocks,const std::array<int,7>&need){
    for(int pairs=0;2*pairs<=blocks&&pairs<=2;++pairs){
        int fixed=blocks-2*pairs;
        if(fixed<0||fixed>5)continue;
        if(pairs==0){if(fixed_sum_codes[fixed].count(encode(need)))return true;}
        else if(pairs==1){
            for(auto code:fixed_sum_codes[fixed])if(pair_unit_codes.count(encode(subv(need,decode(code)))))return true;
        }else{
            for(auto code:fixed_sum_codes[fixed])if(pair_pair_codes.count(encode(subv(need,decode(code)))))return true;
        }
    }
    if(blocks>=4){
        for(auto code:fixed_sum_codes[blocks-4])if(pasch_codes.count(encode(subv(need,decode(code)))))return true;
    }
    return false;
}

static void dfs(int depth){
    if(found)return;++nodes;
    if(!seen[depth].insert(key()).second)return;
    int rem=6-depth,l1=0,coord=-1,mx=0;
    for(int i=0;i<C;++i){int a=std::abs(bal[i]);l1+=a;if(a>mx){mx=a;coord=i;}}
    if(l1>3*rem||mx>rem)return;
    if(rem==1||rem==2){
        std::array<int,7> twice_need{};
        for(int j=0;j<7;++j)twice_need[j]=2*((-total[j]+P)%P)%P;
        const auto code=encode(twice_need);
        if((rem==1&&!pair_unit_codes.count(code))||(rem==2&&!pair_pair_codes.count(code)))return;
    }
    if(coord<0){
        std::array<int,7> need{};for(int j=0;j<7;++j)need[j]=(-total[j]+P)%P;
        if(central_completion(rem,need)){found=true;found_ae=chosen;found_remaining=rem;found_central_need=need;}
        return;
    }
    if(depth==6)return;
    int needed=bal[coord]>0?-1:1;
    for(int id:cand[coord][needed>0?1:0]){
        chosen.push_back(id);for(auto [c,s]:bd[id])bal[c]+=s;for(int j=0;j<7;++j)total[j]=(total[j]+mv[id][j])%P;
        dfs(depth+1);
        for(int j=0;j<7;++j)total[j]=(total[j]-mv[id][j]+P)%P;for(auto [c,s]:bd[id])bal[c]-=s;chosen.pop_back();
        if(found)return;
    }
}

int main(int argc,char**argv){
    std::array<std::array<int,P>,P> ec,es;for(auto&r:ec)r.fill(-2);for(auto&r:es)r.fill(0);int nc=0;
    for(int a=0;a<P;++a)for(int b=a+1;b<P;++b){int na=(-a+P)%P,nb=(-b+P)%P;if(na>nb)std::swap(na,nb);if(a==na&&b==nb){ec[a][b]=ec[b][a]=-1;}else if(ec[a][b]==-2){ec[a][b]=ec[b][a]=nc;es[a][b]=es[b][a]=1;ec[na][nb]=ec[nb][na]=nc;es[na][nb]=es[nb][na]=-1;++nc;}}
    if(nc!=C)return 2;
    std::array<int,P*P*P> idx;idx.fill(-1);
    for(int a=0;a<P;++a)for(int b=a+1;b<P;++b)for(int c=b+1;c<P;++c){int id=tris.size();tris.push_back({a,b,c});idx[(a*P+b)*P+c]=id;}
    bd.resize(tris.size());mv.resize(tris.size());
    for(int id=0;id<(int)tris.size();++id){auto t=tris[id];std::array<std::pair<int,int>,3> ee{{{t[0],t[1]},{t[0],t[2]},{t[1],t[2]}}};std::array<int,C> net{};for(auto[a,b]:ee)if(ec[a][b]>=0)net[ec[a][b]]+=es[a][b];for(int c=0;c<C;++c)if(net[c]){bd[id].push_back({c,net[c]});cand[c][net[c]>0?1:0].push_back(id);}mv[id]=ae_moment(t);if(bd[id].empty())invariant_ids.push_back(id);}
    for(int id=0;id<(int)tris.size();++id){Tri nt=norm_tri(-tris[id][0],-tris[id][1],-tris[id][2]);int ni=idx[(nt[0]*P+nt[1])*P+nt[2]];if(id<ni)nonfixed_reps.push_back(id);}
    fixed_sum_codes[0].insert(0);
    for(int count=1;count<=5;++count)for(auto code:fixed_sum_codes[count-1])for(int id:invariant_ids)fixed_sum_codes[count].insert(encode(addv(decode(code),mv[id])));
    for(int id=0;id<(int)tris.size();++id){std::array<int,7>u{};for(int j=0;j<7;++j)u[j]=2*mv[id][j]%P;auto code=encode(u);if(pair_unit_codes.insert(code).second)pair_unit_vectors.push_back(u);}
    for(int i=0;i<(int)pair_unit_vectors.size();++i)for(int j=i;j<(int)pair_unit_vectors.size();++j)pair_pair_codes.insert(encode(addv(pair_unit_vectors[i],pair_unit_vectors[j])));
    const std::array<std::array<int,6>,7> pis{{
        {{1,0,3,2,5,4}},{{2,4,0,5,1,3}},{{3,4,5,0,1,2}},{{4,5,3,2,0,1}},
        {{5,2,1,4,3,0}},{{5,3,4,1,2,0}},{{5,4,3,2,1,0}}
    }};
    const std::array<Tri,4> pb{{{{0,1,2}},{{0,3,4}},{{1,3,5}},{{2,4,5}}}};
    for(const auto&pi:pis){
        std::vector<std::pair<int,int>> cyc;for(int i=0;i<6;++i)if(i<pi[i])cyc.push_back({i,pi[i]});
        for(int a=0;a<P;++a)for(int b=0;b<P;++b)for(int c=0;c<P;++c){std::array<int,3>par{{a,b,c}};std::array<int,6>lab{};for(int z=0;z<3;++z){lab[cyc[z].first]=par[z];lab[cyc[z].second]=(-par[z]+P)%P;}std::array<int,7>sum{};bool ok=true;for(auto block:pb){Tri u=norm_tri(lab[block[0]],lab[block[1]],lab[block[2]]);if(u[0]==u[1]||u[1]==u[2]){ok=false;break;}int id=idx[(u[0]*P+u[1])*P+u[2]];sum=addv(sum,mv[id]);}if(ok)pasch_codes.insert(encode(sum));}
    }
    std::set<Tri> reps;
    for(int a=0;a<P;++a)for(int b=a+1;b<P;++b)for(int c=0;c<P;++c)if(c!=a&&c!=b){Tri best{a,b,c};for(int s=1;s<P;++s){int aa=a*s%P,bb=b*s%P;if(aa>bb)std::swap(aa,bb);Tri u{aa,bb,c*s%P};if(u<best)best=u;}reps.insert(best);}
    std::vector<Tri> repv(reps.begin(),reps.end());
    int begin=argc>1?std::stoi(argv[1]):0,end=argc>2?std::stoi(argv[2]):(int)repv.size();
    if(begin<0||end<begin||end>(int)repv.size())return 3;
    std::cerr<<"compact_orbits="<<reps.size()<<" range="<<begin<<":"<<end<<" invariant="<<invariant_ids.size()<<" fixed_sums=";for(int i=0;i<=5;++i)std::cerr<<fixed_sum_codes[i].size()<<",";std::cerr<<" pair_units="<<pair_unit_vectors.size()<<" pair_pairs="<<pair_pair_codes.size()<<" pasch="<<pasch_codes.size()<<"\n";
    for(int compact_index=begin;compact_index<end;++compact_index){
        const auto&cp=repv[compact_index];
        ++compact_done;for(auto&x:seen)x.clear();bal.fill(0);chosen.clear();total=compact_moment_vec(cp);
        std::array<std::tuple<int,int,int>,3> ce{{{cp[0],cp[1],1},{cp[0],cp[2],-1},{cp[1],cp[2],-1}}};
        for(auto[a,b,w]:ce)if(ec[a][b]>=0)bal[ec[a][b]]+=w*es[a][b];
        bool zero=true;for(int x:bal)zero&=(x==0);if(zero){std::cout<<"CENTERED_CERTIFIED_INDEX="<<compact_index<<" compact="<<cp[0]<<","<<cp[1]<<","<<cp[2]<<"\n";continue;}
        dfs(0);
        if(found){found_compact=cp;std::cout<<"SAT_INDEX="<<compact_index<<"\n";break;}
        std::cout<<"UNSAT_INDEX="<<compact_index<<" compact="<<cp[0]<<","<<cp[1]<<","<<cp[2]<<" cumulative_nodes="<<nodes<<"\n";
    }
    std::cout<<"SUMMARY range="<<begin<<":"<<end<<" compact_done="<<compact_done<<" nodes="<<nodes<<" found="<<found<<"\n";
    if(found){std::cout<<"compact="<<found_compact[0]<<","<<found_compact[1]<<","<<found_compact[2]<<"\nDEFECT_AE=";for(int id:found_ae)std::cout<<tris[id][0]<<","<<tris[id][1]<<","<<tris[id][2]<<";";std::cout<<" remaining="<<found_remaining<<" need=";for(int x:found_central_need)std::cout<<x<<",";std::cout<<"\n";}
}
