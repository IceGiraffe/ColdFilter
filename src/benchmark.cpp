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
#include <fstream>

#include <sys/stat.h>

#include "./On-Off-Sketch/PE/OO_PE.h"
#include "./On-Off-Sketch/FPI/OO_FPI.h"
#include "SF_noSIMD.h"

using namespace std;

using COUNT_TYPE = uint32_t;
using DATA_TYPE = uint32_t;
typedef std::unordered_map<DATA_TYPE, COUNT_TYPE> HashMap;

HashMap mp;
vector<uint32_t> insert_data;

int TOTAL = 0;
int T = 0;
const int _T = 1600;
int LENGTH;
vector<HashMap *> WindowCount;

void load_file(const char * PATH){
    FILE* file = fopen(PATH, "rb");
    COUNT_TYPE number = 0;
    DATA_TYPE item;
    HashMap record;
    

    TOTAL = 0;
    T = 0;

    struct stat statbuf;
    stat(PATH, &statbuf);
    LENGTH = ceil(statbuf.st_size / (double)(_T * (13 + 8)));
    
    char tmp[13];
    char ts[8];
    while(fread(&tmp, 13, 1, file) > 0){
        item = *(uint32_t *)tmp;
        int rsize = fread(ts, 1, 8, file);
        if(rsize != 8) {
            std::cout << "TimeStamp Failed.." << std::endl;
            break;
        }
        
        if(number % LENGTH == 0){
            // auto h = new HashMap;
            // WindowCount.push_back(h);
            T += 1;
            // std::cout << "Period: " << T << " Begin..." << std::endl;
        }
        number += 1;
        insert_data.push_back(item);
        // (*WindowCount[T-1])[item] += 1;

        if(record[item] != T){
            record[item] = T;
            mp[item] += 1;
            TOTAL += 1;
        }
    }

    cout << "Total Number: " << insert_data.size() << endl;
    cout << "Total Window: " << T << endl;  
    cout << "Window Length: " << LENGTH << endl;
    cout << "Distinct Item: " << mp.size() << endl;
    uint32_t max_freq = 0;
    for (auto itr: mp) {
        max_freq = max(max_freq, itr.second);
    }
    printf("Max frequency : %d\n", max_freq);
    fclose(file);
}

void check_dis(){
    fstream fout("persistent.txt");
    for(auto it = mp.begin();it != mp.end();++it){
        if(it->second > 100){
            int sum = 0;
            for(int i = 0; i < T; i++){
                if((*WindowCount[i])[it->first] == 1){
                    sum ++;
                }
            }
            fout << it->first << " " << it->second << " " << sum << endl;
        }
        // number[pos] += 1;
    }
}

void PE(){
    int section = 10;
    // 500 KB
    cout << "====================PureOO_PE=================" << endl;
    auto MyOO = new OO_PE<DATA_TYPE, COUNT_TYPE> (2, 500000 / 2.0 / (BITSIZE + sizeof(COUNT_TYPE)));
    
    // insert items to MyOO
    {
        DATA_TYPE item;
        COUNT_TYPE number = 0, windowId = 0;

        for(auto item: insert_data){
            if(number % LENGTH == 0){
                windowId += 1;
                MyOO->NewWindow(windowId);
            }
            number += 1;

            MyOO->Insert(item, windowId);
        }
    }
    // check the result
    {
        uint32_t* aae = new uint32_t[section];
        uint32_t* number = new uint32_t[section];
        uint32_t total_aae = 0;
        
        memset(aae, 0, sizeof(uint32_t) * section);
        memset(number, 0, sizeof(uint32_t) * section);
        cout << "T " << T << endl;
        for(auto it = mp.begin();it != mp.end();++it){
            COUNT_TYPE value = MyOO->Query(it->first);
            uint32_t pos = (it->second * section - 1) / T;
            if(abs(int32_t(it->second - value)) > 1600) cout << it->second << " " << it->first << " " << value << endl;
            total_aae += abs(int32_t(it->second - value));
            aae[pos] += abs(int32_t(it->second - value));
            number[pos] += 1;
        }

        for(uint32_t i = 0;i < section;++i){
            // printf("Part: %d with %d Item\n", i, number[i]);
            if(number[i] != 0)
                std::cout << aae[i] / (double)number[i] << endl;
            else
                std::cout << "NULL" << endl;
        }

        delete [] aae;
        delete [] number;
        std::cout << "AAE: " << total_aae / mp.size() << std::endl;

    }
    
    cout << "====================SF+OO_PE==================" << endl;
    auto Sf = new SlidingFilter(64, 0, 256, 4, 2, 10000, 10000, 190, 20);
    auto SfOO = new OO_PE<DATA_TYPE, COUNT_TYPE> (2, 420000 / 2.0 / (BITSIZE + sizeof(COUNT_TYPE)));
    // insert items to Sf
    {
        DATA_TYPE item;
        COUNT_TYPE number = 0, windowId = 0;
        int SF_Number = 0;
        bool flag = 1;
        for(auto item: insert_data){
            if(number % LENGTH == 0){
                flag = 1;
                // for(int i = 0; i < Sf->bucket_num1 ; i++){
                    // for(int j = 0; j < Sf->cols; j++){
                        // if(Sf->Bucket1[i].fp[j] != 0)
                        //     SfOO->Insert(Sf->Bucket1[i].fp[j], windowId);
                        // Sf->Bucket1[i].fp[j] = 0;
                        // Sf->Bucket1[i].counter[j] = 0;
                //     }
                // }
                windowId += 1;
                SfOO->NewWindow(windowId);
            }
            number += 1;

            int res = Sf->insert(item, 1);
            // if(item == 3918694019 && flag) {
            //     flag = 0;
            //     cout << res << endl;
            // }
            if(res >= 10){
                SF_Number++;
                SfOO->Insert(item, windowId);
            }
        }
        cout << "SF_Number: "<< SF_Number << endl;
    }
    // check the result
    {
        uint32_t* aae = new uint32_t[section];
        uint32_t* number = new uint32_t[section];
        uint32_t total_aae = 0;
        
        memset(aae, 0, sizeof(uint32_t) * section);
        memset(number, 0, sizeof(uint32_t) * section);
        cout << "T " << T << endl;
        for(auto it = mp.begin();it != mp.end();++it){
            COUNT_TYPE value = SfOO->Query(it->first);
            uint32_t pos = (it->second * section - 1) / T;
            // if(abs(int32_t(it->second - value)) > 100) cout << it->second << " " << it->first << " " << value << endl;
            total_aae += abs(int32_t(it->second - value));
            aae[pos] += abs(int32_t(it->second - value));
            number[pos] += 1;
        }

        for(uint32_t i = 0;i < section;++i){
            // printf("Part: %d with %d Item\n", i, number[i]);
            if(number[i] != 0)
                std::cout << aae[i] / (double)number[i] << endl;
            else
                std::cout << "NULL" << endl;
        }

        delete [] aae;
        delete [] number;
        std::cout << "AAE: " << total_aae / mp.size() << std::endl;
    }
}

void FPI(){
    // 200 KB
    int HIT = TOTAL * 0.00005;
    cout << "HIT: " << HIT << endl;
    cout << "====================PureOO_FPI=================" << endl;
    auto MyOO = new OO_FPI<DATA_TYPE, COUNT_TYPE, 8> (150000);
    cout << MyOO->length * 8<< endl;
    // insert items to MyOO
    {
        DATA_TYPE item;
        COUNT_TYPE number = 0, windowId = 0;

        for(auto item: insert_data){
            if(number % LENGTH == 0){
                windowId += 1;
                MyOO->NewWindow(windowId);
            }
            number += 1;

            MyOO->Insert(item, windowId);
        }
    }
    // check the result
    {
        double real = 0, estimate = 0, both = 0;
        double aae = 0, cr = 0, pr = 0, f1 = 0;

        for(auto it = mp.begin();it != mp.end();++it){
            COUNT_TYPE value = MyOO->Query(it->first);

            if(value > HIT){
                estimate += 1;
                if(it->second > HIT) {
                    both += 1;
                    aae += abs(int32_t(it->second) - int32_t(value));
                }
            }
            if(it->second > HIT)
                real += 1;
        }

        if(both <= 0){
            std::cout << "Not Find Real Persistent" << std::endl;
        }
        else{
            aae /= both;
        }

        cr = both / real;

        if(estimate <= 0){
            std::cout << "Not Find Persistent" << std::endl;
        }
        else{
            pr = both / estimate;
        }

        if(cr == 0 && pr == 0)
            f1 = 0;
        else
            f1 = (2*cr*pr)/(cr+pr);

	    cout << HIT << endl;
        printf("Estimate: %f, Both: %f, Real: %f\n", estimate, both, real);
        printf("%f %f %f %f\n", aae, pr, cr, f1);
        cout << "AAE: " << aae << endl;
        cout << cr << " " << pr << endl;
        cout << "F1: " << f1 << endl;
    }
    
    cout << "====================SF+OO_FPI==================" << endl;
    auto Sf = new SlidingFilter(128, 0, 16, 4, 2, 10000, 10000, 190, 20);
    auto SfOO = new OO_FPI<DATA_TYPE, COUNT_TYPE, 8> (100000 - 128 * 16 * 6);
    // insert items to Sf
    {
        DATA_TYPE item;
        COUNT_TYPE number = 0, windowId = 0;
        int SF_Number = 0;
        bool flag = 1;
        for(auto item: insert_data){
            if(number % LENGTH == 0){
                flag = 1;
                for(int i = 0; i < Sf->bucket_num1 ; i++){
                    for(int j = 0; j < Sf->cols; j++){
                        // if(Sf->Bucket1[i].fp[j] != 0)
                            // SfOO->Insert(Sf->Bucket1[i].fp[j], windowId);
                        // Sf->Bucket1[i].fp[j] = 0;
                        // Sf->Bucket1[i].counter[j] = 0;
                    }
                }
                windowId += 1;
                SfOO->NewWindow(windowId);
            }
            number += 1;

            int res = Sf->insert(item, 1);
            // if(item == 3918694019 && flag) {
            //     flag = 0;
            //     cout << res << endl;
            // }
            if(res >= 2){
                SF_Number++;
                SfOO->Insert(item, windowId);
                
            }
        }
        cout << "SF_Number: "<< SF_Number << endl;
    }
    // check the result
    {
        double real = 0, estimate = 0, both = 0;
        double aae = 0, cr = 0, pr = 0, f1 = 0;

        for(auto it = mp.begin();it != mp.end();++it){
            COUNT_TYPE value = SfOO->Query(it->first);

            if(value > HIT){
                estimate += 1;
                if(it->second > HIT) {
                    both += 1;
                    aae += abs(int32_t(it->second) - int32_t(value));
                }
            }
            if(it->second > HIT)
                real += 1;
        }

        if(both <= 0){
            std::cout << "Not Find Real Persistent" << std::endl;
        }
        else{
            aae /= both;
        }

        cr = both / real;

        if(estimate <= 0){
            std::cout << "Not Find Persistent" << std::endl;
        }
        else{
            pr = both / estimate;
        }

        if(cr == 0 && pr == 0)
            f1 = 0;
        else
            f1 = (2*cr*pr)/(cr+pr);

	    cout << HIT << endl;
        printf("Estimate: %f, Both: %f, Real: %f\n", estimate, both, real);
        printf("%f %f %f %f\n", aae, pr, cr, f1);
        cout << "AAE: " << aae << endl;
        cout << cr << " " << pr << endl;
        cout << "F1: " << f1 << endl;
    }
}

int main(){
    load_file("/root/CAIDA/130000.dat");
    // check_dis();
    // PE();
    FPI();
    return 0;
}