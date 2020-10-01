/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MyUtils.cpp
 * Author: sylwester
 * 
 * Created on November 28, 2018, 1:59 PM
 */

#include <Utils/MyUtils.h>
#include <Params.h>
#include <Global.h>

#include "Utils/MyUtils.h"

//MyUtils::MyUtils() {
//}
//
//MyUtils::MyUtils(const MyUtils& orig) {
//}
//
//MyUtils::~MyUtils() {
//}

int MyUtils::parseInt(string s) {
    stringstream ss(s);
    int a;
    ss >> a;
    return a;
}

string MyUtils::getComplimentaryString(string s) {
    string res = s;
    for (int i = 0; i < s.size(); i++) {
        if (Params::getNukl(s[i]) == Params::A) res[i] = 'T';
        if (Params::getNukl(s[i]) == Params::C) res[i] = 'G';
        if (Params::getNukl(s[i]) == Params::G) res[i] = 'C';
        if (Params::getNukl(s[i]) == Params::T) res[i] = 'A';
    }
//    string res = "";
//    for( int i=0; i<s.size(); i++ ){
//        if( Params::getNukl(s[i]) == Params::A ) res += Params::getNuklAsString( Params::T );
//        if( Params::getNukl(s[i]) == Params::C ) res += Params::getNuklAsString( Params::G );
//        if( Params::getNukl(s[i]) == Params::T ) res += Params::getNuklAsString( Params::A );
//        if( Params::getNukl(s[i]) == Params::G ) res += Params::getNuklAsString( Params::C );
//    }
    return res;

}

int MyUtils::getNearestLowerPrime(int N) {
    int X = N - 1;
    bool prime = false;
    while (!prime) {
        prime = true;

        if (X % 2 == 0) {
            prime = false;
            X--;
        } else {
            for (int i = 3; i * i <= X; i += 2) {
                if (X % i == 0) {
                    prime = false;
                    break;
                }
            }

            if (!prime) X--;
            else break;
        }

    }
    return X;
}


void MyUtils::process_mem_usage(double vm_usage, double resident_set) {
    vm_usage = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
            >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
            >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;

    cerr << endl;
    DEBUG(vm_usage);
    DEBUG(resident_set);
    cerr << endl;
}