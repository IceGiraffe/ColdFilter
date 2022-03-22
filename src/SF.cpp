#include <iostream>
#include <string>
#include <queue>
#include <stack>
#include <list>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <memory>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <cstring>

#include "SpaceSaving.h"
#include "SC_SpaceSaving.h"

// Bob hash
#include "BOBHash32.h"
#include "CUSketch.h"
#include "CUHeap.h"
#include "SPA.h"
// SIMD 
#include <immintrin.h>
#define UNIX
#ifdef UNIX
#include <x86intrin.h>
#else
#include <intrin.h>
#endif

using namespace std;

#define MAX_INSERT_PACKAGE 32000000

unordered_map<uint32_t, int> ground_truth;
uint32_t insert_data[MAX_INSERT_PACKAGE];
uint32_t query_data[MAX_INSERT_PACKAGE];

template<uint32_t waitingSize>
class WaitingRoom{
public:
    uint32_t fingerPrint[waitingSize] __attribute__ ((aligned (16)));
    uint32_t counter[waitingSize] __attribute__ ((aligned (16)));
    BOBHash32 * hasher;
    WaitingRoom(){
        memset(fingerPrint, 0, sizeof(fingerPrint));
        memset(counter, 0, sizeof(counter));
        hasher = new BOBHash32(18973);
    }

    void insert(uint32_t key, uint32_t f){
        auto fp = hasher->run((const char *)&key, 4);
        auto fp_packed = _mm_set1_epi32((int)fp);
        int matched = 0;
        if (waitingSize == 16) {
            __m128i *keys_p = (__m128i *)fingerPrint;

            __m128i a_comp = _mm_cmpeq_epi32(fp_packed, keys_p[0]);
            __m128i b_comp = _mm_cmpeq_epi32(fp_packed, keys_p[1]);
            __m128i c_comp = _mm_cmpeq_epi32(fp_packed, keys_p[2]);
            __m128i d_comp = _mm_cmpeq_epi32(fp_packed, keys_p[3]);

            a_comp = _mm_packs_epi32(a_comp, b_comp);
            c_comp = _mm_packs_epi32(c_comp, d_comp);
            a_comp = _mm_packs_epi32(a_comp, c_comp);

            matched = _mm_movemask_epi8(a_comp);
        } else if (CounterNum == 4) {
            __m128i *keys_p = (__m128i *)fingerPrint;
            __m128i a_comp = _mm_cmpeq_epi32(fp_packed, keys_p[0]);
            matched = _mm_movemask_ps(*(__m128 *)&a_comp);
        } else {
            throw std::logic_error("Not implemented.");
        }

        if(matched){
            int matched_index = _tzcnt_u32((uint32_t) matched);
            counter[matched_index] += f;
            return true;
        }
    }

    void check(uint32_t key){
        auto fp = hasher->run((const char *)&key, 4);
        auto fp_packed = _mm_set1_epi32((int)fp);
        int matched = 0;
        if (waitingSize == 16) {
            __m128i *keys_p = (__m128i *)fingerPrint;

            __m128i a_comp = _mm_cmpeq_epi32(fp_packed, keys_p[0]);
            __m128i b_comp = _mm_cmpeq_epi32(fp_packed, keys_p[1]);
            __m128i c_comp = _mm_cmpeq_epi32(fp_packed, keys_p[2]);
            __m128i d_comp = _mm_cmpeq_epi32(fp_packed, keys_p[3]);

            a_comp = _mm_packs_epi32(a_comp, b_comp);
            c_comp = _mm_packs_epi32(c_comp, d_comp);
            a_comp = _mm_packs_epi32(a_comp, c_comp);

            matched = _mm_movemask_epi8(a_comp);
        } else if (CounterNum == 4) {
            __m128i *keys_p = (__m128i *)fingerPrint;
            __m128i a_comp = _mm_cmpeq_epi32(fp_packed, keys_p[0]);
            matched = _mm_movemask_ps(*(__m128 *)&a_comp);
        } else {
            throw std::logic_error("Not implemented.");
        }
        if(matched){
            int matched_index = _tzcnt_u32((uint32_t) matched);
            fingerPrint[matched_index] = 0;
            return counter[matched_index];
        }
    }
};


// multi-layer to avoid full and record aged hot items
template<uint32_t L1BucketNum, uint32_t L2BucketNum, uint32_t CounterNum, uint32_t L1_THRES, uint32_t ss_capacity>
class SlidingFilter{
public:
    // CUSketch<int(1024 * 5), 3> ss;
    SpaceSaving<ss_capacity> ss;
    uint32_t curTick;
    BOBHash32 * l1Hash;
    BOBHash32 * l2Hash;

    // Number of Insertions to Space Saving
    int exceed_cnt = 0;

    void printBasicInfo(){
        int total_bytes = 4 * 3 * CounterNum * (L1BucketNum + L2BucketNum);
        printf("Sliding Filter Memory Usage:\n");
        printf("\tTotal memory usage is %d bytes = %f KB.\n", total_bytes, total_bytes / 1024.0);
        printf("\tSpace Saving Capacity is %d\n", ss_capacity);
    }

    struct Bucket
    {
        /* data */
        uint32_t fingerPrint[CounterNum] __attribute__ ((aligned (16)));
        uint32_t timeStamp[CounterNum] __attribute__ ((aligned (16)));
        uint32_t counter[CounterNum];
        uint32_t pos;
        Bucket(): pos(0){
            memset(fingerPrint, 0, sizeof(fingerPrint));
            memset(timeStamp, 0, sizeof(timeStamp));
            memset(counter, 0, sizeof(counter));
        }
    };

    Bucket L1Cache[L1BucketNum] __attribute__ ((aligned (16)));

    // Bucket L2Cache[L2BucketNum] __attribute__ ((aligned (16)));

    // Constructor
    SlidingFilter():curTick(0){
        l1Hash = new BOBHash32(CounterNum + 1);
        l2Hash = new BOBHash32(CounterNum + 2);
    }

    // Destructor
    ~SlidingFilter(){
        delete l1Hash;
        delete l2Hash;
    }

    // get timeStamp
    inline uint32_t tick(){
        return ++curTick;
    }

    // quick modular
    inline uint32_t multiply_high_u32(uint32_t x, uint32_t y) const {
        return (uint32_t) (((uint64_t) x * (uint64_t) y) >> 32);
    }

    // insert a key
    void insert(uint32_t key, uint32_t f = 1){
        bool succ = insertL1(key, f);
        return;
        // if(succ) return;
        // insertL2(key, f);
    }

    inline void build(uint32_t * items, int n)
    {
        for (int i = 0; i < n; ++i) {
            this->insert(items[i], 1);
        }
        this->synch();
    }

    void get_top_k(uint32_t k, uint32_t items[])
    {
        ss.get_top_k(k, items);
    }

    void get_top_k_with_frequency(uint32_t k, vector<pair<uint32_t, uint32_t>> & result)
    {
        ss.get_top_k_with_frequency(k, result);
    }

    // insert a key to L1
    bool insertL1(uint32_t key, uint32_t f){
        uint32_t ts = tick();
        auto key_packed = _mm_set1_epi32((int)key);

        int L1Pos = multiply_high_u32( l1Hash->run((const char *)&key, 4), L1BucketNum);
        Bucket & hashed_bucket = L1Cache[L1Pos];

        int matched = 0;
        if (CounterNum == 16) {
            __m128i *keys_p = (__m128i *)hashed_bucket.fingerPrint;

            __m128i a_comp = _mm_cmpeq_epi32(key_packed, keys_p[0]);
            __m128i b_comp = _mm_cmpeq_epi32(key_packed, keys_p[1]);
            __m128i c_comp = _mm_cmpeq_epi32(key_packed, keys_p[2]);
            __m128i d_comp = _mm_cmpeq_epi32(key_packed, keys_p[3]);

            a_comp = _mm_packs_epi32(a_comp, b_comp);
            c_comp = _mm_packs_epi32(c_comp, d_comp);
            a_comp = _mm_packs_epi32(a_comp, c_comp);

            matched = _mm_movemask_epi8(a_comp);
        } else if (CounterNum == 4) {
            __m128i *keys_p = (__m128i *)hashed_bucket.fingerPrint;
            __m128i a_comp = _mm_cmpeq_epi32(key_packed, keys_p[0]);
            matched = _mm_movemask_ps(*(__m128 *)&a_comp);
        } else {
            throw std::logic_error("Not implemented.");
        }
        
        // match L1 Cache
        if(matched){
            int matched_index = _tzcnt_u32((uint32_t) matched);
            hashed_bucket.counter[matched_index] += f;
            hashed_bucket.timeStamp[matched_index] = ts;
            return true;
        }

        if(hashed_bucket.pos != CounterNum){
            int pos = hashed_bucket.pos;
            hashed_bucket.fingerPrint[pos] = key;
            hashed_bucket.counter[pos] = f;
            hashed_bucket.timeStamp[pos] = ts;
            hashed_bucket.pos ++;
            return true;
        }

        // find the agest one
        int pos = 0;
        int minV = hashed_bucket.timeStamp[0];
        for(int i = 1; i < CounterNum; i++){
            if(hashed_bucket.timeStamp[i] < minV){
                minV = hashed_bucket.timeStamp[i];
                pos = i;
            }
        }
        auto old_key = hashed_bucket.fingerPrint[pos];
        auto old_c = hashed_bucket.counter[pos];

        hashed_bucket.fingerPrint[pos] = key;
        hashed_bucket.counter[pos] = f;
        hashed_bucket.timeStamp[pos] = ts;
        
        if(old_c >= L1_THRES) {
            ++exceed_cnt;
            // cout << "Exceed " << exceed_cnt  << endl;
            ss.insert(old_key, old_c);
        }
        return false;
    }

    /* query a key from CU
    uint32_t query(uint32_t key){
        return 1;
    }
    */

    void synch(){
        for(int i = 0; i < L1BucketNum; i++){
            for(int j = 0; j < CounterNum; j++){
                if(L1Cache[i].counter[j] >= L1_THRES)
                    ss.insert(L1Cache[i].fingerPrint[j], L1Cache[i].counter[j]);
            }
        }
    }


    /* Query a key in L1
    // query a key in L1
    uint32_t queryL1(uint32_t key){
        uint32_t ts = tick();
        auto key_packed = _mm_set1_epi32(key);

        int L1Pos = multiply_high_u32( l1Hash->run((const char *)&key, 4), L1BucketNum);
       
        // cout << L1Pos << endl;

        Bucket & hashed_bucket = L1Cache[L1Pos];
        int matched = 0;
        if (CounterNum == 16) {
            __m128i *keys_p = (__m128i *)hashed_bucket.fingerPrint;

            __m128i a_comp = _mm_cmpeq_epi32(key_packed, keys_p[0]);
            __m128i b_comp = _mm_cmpeq_epi32(key_packed, keys_p[1]);
            __m128i c_comp = _mm_cmpeq_epi32(key_packed, keys_p[2]);
            __m128i d_comp = _mm_cmpeq_epi32(key_packed, keys_p[3]);

            a_comp = _mm_packs_epi32(a_comp, b_comp);
            c_comp = _mm_packs_epi32(c_comp, d_comp);
            a_comp = _mm_packs_epi32(a_comp, c_comp);

            matched = _mm_movemask_epi8(a_comp);
        } else if (CounterNum == 4) {
            __m128i *keys_p = (__m128i *)hashed_bucket.fingerPrint;
            __m128i a_comp = _mm_cmpeq_epi32(key_packed, keys_p[0]);
            matched = _mm_movemask_ps(*(__m128 *)&a_comp);
        } else {
            throw std::logic_error("Not implemented.");
        }
        
        // match L1 Cache
        if(matched){
            int matched_index = _tzcnt_u32((uint32_t) matched);
            hashed_bucket.timeStamp[matched_index] = ts;
            return hashed_bucket.counter[matched_index];
        }
        else return 0;
    }
    */


    /* insert a key to L2
    // insert a key to L2
    void insertL2(uint32_t key, uint32_t f){
        uint32_t ts = tick();
        auto key_packed = _mm_set1_epi32(key);

        int L2Pos = multiply_high_u32( l2Hash->run((const char *)&key, 4), L2BucketNum);

        // cout << L2Pos << endl;

        Bucket & hashed_bucket = L2Cache[L2Pos];
        assert(CounterNum == 4);
        __m128i *fp_packed = (__m128i *)hashed_bucket.fingerPrint;
        __m128i comp = _mm_cmpeq_epi32(key_packed, fp_packed[0]);
        auto matched = _mm_movemask_ps(*(__m128 *)&comp);
        
        // match L2 Cache
        if(matched){
            int matched_index = _tzcnt_u32((uint32_t) matched);
            hashed_bucket.counter[matched_index] += f;
            hashed_bucket.timeStamp[matched_index] = ts;
            return;
        }

        if(hashed_bucket.pos != CounterNum){
            int pos = hashed_bucket.pos;
            hashed_bucket.fingerPrint[pos] = key;
            hashed_bucket.counter[pos] = f;
            hashed_bucket.timeStamp[pos] = ts;
            hashed_bucket.pos ++;
            return;
        }

        // find the agest one
        int pos = 0;
        int minV = hashed_bucket.timeStamp[0];
        for(int i = 1; i < CounterNum; i++){
            if(hashed_bucket.timeStamp[i] < minV){
                minV = hashed_bucket.timeStamp[i];
                pos = i;
            }
        }
        
        auto old_key = hashed_bucket.fingerPrint[pos];
        auto old_c = hashed_bucket.counter[pos];
        ss.insert(old_key, old_c);
        
        hashed_bucket.fingerPrint[pos] = key;
        hashed_bucket.counter[pos] = f;
        hashed_bucket.timeStamp[pos] = ts;
        
        return;
    }
    */


    /* Query L2
    // query a key in L2 
    uint32_t queryL2(uint32_t key){
        // cout << "Query L2..." << endl;
        uint32_t ts = tick();
        auto key_packed = _mm_set1_epi32(key);

        int L2Pos = multiply_high_u32( l2Hash->run((const char *)&key, 4), L2BucketNum);
       
        // cout << L2Pos << endl;

        Bucket & hashed_bucket = L2Cache[L2Pos];
        assert(CounterNum == 4);
        __m128i *fp_packed = (__m128i *)hashed_bucket.fingerPrint;
        __m128i comp = _mm_cmpeq_epi32(key_packed, fp_packed[0]);
        auto matched = _mm_movemask_ps(*(__m128 *)&comp);
        
        // match L2 Cache
        if(matched){
            int matched_index = _tzcnt_u32((uint32_t) matched);
            hashed_bucket.timeStamp[matched_index] = ts;
            return hashed_bucket.counter[matched_index];
        }
        else return 0;
    }
    */
};

// load data
int load_data(const char *filename) {
    FILE *pf = fopen(filename, "rb");
    if (!pf) {
        cerr << filename << " not found." << endl;
        exit(-1);
    }

    ground_truth.clear();

    char ip[13];
    char ts[8];
    int ret = 0;
    while (1) {
        size_t rsize;
        rsize = fread(ip, 1, 13, pf);
        if(rsize != 13) break;
        rsize = fread(ts, 1, 8, pf);
        if(rsize != 8) break;

        uint32_t key = *(uint32_t *) ip;
        insert_data[ret] = key;
        ground_truth[key]++;
        ret++;
        if (ret == MAX_INSERT_PACKAGE){
            cout << "MAX_INSERT_PACKAGE" << endl;
            break;
        }
    }
    fclose(pf);

    int i = 0;
    for (auto itr: ground_truth) {
        query_data[i++] = itr.first;
    }

    printf("Total stream size = %d\n", ret);
    printf("Distinct item number = %ld\n", ground_truth.size());

    int max_freq = 0;
    for (auto itr: ground_truth) {
        max_freq = std::max(max_freq, itr.second);
    }
    printf("Max frequency = %d\n", max_freq);

    return ret;
}


pair<double, double> ss_compare_value_with_ground_truth(uint32_t k, vector<pair<uint32_t, uint32_t>> & result)
{
    // prepare top-k ground truth
    vector<pair<uint32_t, int>> gt_ordered;
    for (auto itr: ground_truth) {
        gt_ordered.emplace_back(itr);
    }
    std::sort(gt_ordered.begin(), gt_ordered.end(), [](const std::pair<uint32_t, int> &left, const std::pair<uint32_t, int> &right) {
        return left.second > right.second;
    });
    set<uint32_t> set_gt;
    int i = 0;
    int th;
    for (auto itr: gt_ordered) {
        if (i >= k && itr.second < th) {
            break;
        }
        set_gt.insert(itr.first);
        i++;
        if (i == k) {
            th = itr.second;
        }
    }

    double aae = 0;

    set<uint32_t> set_rp;
    unordered_map<uint32_t, uint32_t> mp_rp;

    int lp = 0;
    for (lp = 0; lp < k; ++lp) {
        set_rp.insert(result[lp].first);
        mp_rp[result[lp].first] = result[lp].second;
    }

    vector<uint32_t> intersection(k);
    auto end_itr = std::set_intersection(
            set_gt.begin(), set_gt.end(),
            set_rp.begin(), set_rp.end(),
            intersection.begin()
    );

    for (auto itr = intersection.begin(); itr != end_itr; ++itr) {
        int diff = int(mp_rp[*itr]) - int(ground_truth[*itr]);
        aae += abs(diff);
    }

    int num = end_itr - intersection.begin();
    num = num > 0 ? num : 1;

    return make_pair(double(num) / k, aae / num);
}

template<uint32_t topk, uint32_t L2Thres>
void demo_ss(int packet_num)
{
    printf("\nExp for top-k query:\n");
    // Sliding Filter: Memory = 200KB, H = 2.5K = 2560
    // ColdFilter: Memory = 200KB, H = 2.5K = 2560
    // Space Saving: H = 4608
    /*
     * Total stream size = 27121713
     * Distinct item number = 253906
     * Max frequency = 916516
     */
    printf("\ttest top %d\n", topk);

    auto sf_ss = new SlidingFilter<1024, 0, 16, 400, 2560>();
    auto cf_ss = new SC_SpaceSaving<2560, 200*1024, L2Thres, 16>();
    auto ss = new SpaceSaving<4608>();
    
    sf_ss->printBasicInfo();
    cf_ss->sc.print_basic_info();
    
    
    timespec dtime1, dtime2;
    long long delay = 0.0;
    
    int k = topk;
    vector<pair<uint32_t, uint32_t>> result(k);
    
    printf("-----------------------------------------------\n");
    
    vector< pair<double, double> > res;
    // SlidingFilter
    clock_gettime(CLOCK_MONOTONIC, &dtime1);
    sf_ss->build(insert_data, packet_num);
    clock_gettime(CLOCK_MONOTONIC, &dtime2);
    delay = (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
    cout << delay << endl;
    
    {
        sf_ss->get_top_k_with_frequency(k, result);
        auto ret = ss_compare_value_with_ground_truth(k, result);
        res.push_back(ret);
        printf("\tSliding Filter \n");
        printf("\t  accuracy %lf, AAE %lf\n", ret.first, ret.second);
    }
    
    // return;
    // ColdFilter
    clock_gettime(CLOCK_MONOTONIC, &dtime1);
    cf_ss->build(insert_data, packet_num);
    clock_gettime(CLOCK_MONOTONIC, &dtime2);
    delay = (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
    cout << delay << endl;
    
    {
        cf_ss->get_top_k_with_frequency(k, result);
        auto ret = ss_compare_value_with_ground_truth(k, result);
        res.push_back(ret);
        printf("\tCold Filter:\n");
        printf("\t  accuracy %lf, AAE %lf\n", ret.first, ret.second);
    }

    // SpaceSaving
    clock_gettime(CLOCK_MONOTONIC, &dtime1);
    ss->build(insert_data, packet_num);
    clock_gettime(CLOCK_MONOTONIC, &dtime2);
    delay = (long long)(dtime2.tv_sec - dtime1.tv_sec) * 1000000000LL + (dtime2.tv_nsec - dtime1.tv_nsec);
    cout << delay << endl;
    
    {
        ss->get_top_k_with_frequency(k, result);
        auto ret = ss_compare_value_with_ground_truth(k, result);
        res.push_back(ret);
        printf("\tSpaceSaving:\n");
        printf("\t  accuracy %lf, AAE %lf\n", ret.first, ret.second);
    }

    for(auto k: res){
        cout << k.first << endl;
    }
    return;
}

/* Used with SF+CU ... depreciated
void myTestTopK(int packet_num){
    vector<pair<uint32_t, int>> topK;
    for (auto itr: ground_truth){
        topK.emplace_back(itr);
    }
    sort(topK.begin(), topK.end(), [](auto & a, auto & b){return a.second > b.second;});
    
    topK.erase(topK.begin() + 100, topK.end());
    

    // new SF    
    auto sf = new SlidingFilter<32, 0, 16, 5, 0>;

    sf->printBasicInfo();

    for (int i = 0; i < packet_num; ++i) {
        sf->insert(insert_data[i]);
    }
    sf->synch();

    long long tot_ae = 0;
    
    for (auto itr: topK) {
        int report_val = sf->query(itr.first);
        cout << report_val << " " << itr.second << endl;
        tot_ae += abs(report_val - itr.second);
    }
    printf("Report\n");
    printf("\tSF AAE: %lf\n", double(tot_ae) / topK.size());

    auto mycuh = new CUHeap<100, 10*1024, 3>;

    for (int i = 0; i < packet_num; ++i) {
        mycuh->insert(insert_data[i]);
    }

    tot_ae = 0;
    cout << tot_ae << endl;

    vector<std::pair<uint32_t, uint32_t>> cuTopKRes;

    mycuh->get_top_k_with_frequency(100, cuTopKRes);

    for (auto itr: cuTopKRes) {
        int report_val = itr.second;
        cout << report_val << " " << ground_truth[itr.first] << endl;
        tot_ae += abs(ground_truth[itr.first] - report_val);
    }

    printf("Report\n");
    printf("\tHEAP AAE: %lf\n", double(tot_ae) / cuTopKRes.size());
}
*/

int main(){
    // int packet_num = load_data("../data/sample.dat");
    int packet_num = load_data("/root/CAIDA/130000.dat");

    // demo_ss<uint32_t k> k from 32 to 1024
    map<int, int> thres;
    thres[32] = 40237;
    thres[64] = 36593;
    thres[128] = 17856;
    thres[256] = 12028;
    thres[512] = 8279;
    thres[1024] = 4580;
    /*
     * 32 44708 40237
     * 64 29548 26593
     * 128 19841 17856
     * 256 13365 12028
     * 512 9199 8279
     * 1024 5089 4580
     */
    demo_ss<64, 8279>(packet_num);

    return 0;
}