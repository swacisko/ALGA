//
// Created by sylwester on 12/20/18.
//

#ifndef GENOMEALIGNMENT_STATISTICSGENERATORBIGDATA_H
#define GENOMEALIGNMENT_STATISTICSGENERATORBIGDATA_H

#include<iostream>
#include <cmath>
#include "Params.h"
#include<map>

using namespace std;

class StatisticsGeneratorBigData {

public:

    StatisticsGeneratorBigData() {

    }

    template<class _T>
    static void addData(string s, _T val) {
//        if( Params::WRITE_STATISTICS == false ) return;
        if (data.find(s) == data.end()) data[s] = Stat();

        data[s].sum += val;
        data[s].sumOfSquares += (double) val * (double) val;
        data[s].numberOfElements++;
        data[s].maxVal = max(data[s].maxVal, (double) val);
        data[s].minVal = min(data[s].minVal, (double) val);
    }

    static double getAverage(string s) {
        if (data[s].numberOfElements == 0) return 0;
        else
            return data[s].sum / data[s].numberOfElements;
    }

    static double getStandardDeviation(string s) {
        double avg = getAverage(s);

        double res = sqrt(data[s].sumOfSquares / data[s].numberOfElements - avg * avg);
        return res;
    }

    static void writeAllStatistics(string s) {
        cerr << "number of elements:\t" << data[s].numberOfElements << endl;
        cerr << "sum: " << data[s].sum << endl;
        cerr << "sum of squares: " << data[s].sumOfSquares << endl;
        cerr << "average: " << getAverage(s) << endl;
        cerr << "max: " << data[s].maxVal << endl;
        cerr << "min: " << data[s].minVal << endl;
        cerr << "standard deviation: " << getStandardDeviation(s) << endl;
    }

    static void writeAllStatistics() {
        for (auto a : data) {
            cerr << a.first << ":" << endl;
            writeAllStatistics(a.first);
            cerr << endl;
        }
    }

    static string SAME_KMER_SIZES;
    static string SWAT_ALG_BRANCHES;
    static string COMPARISONS_IN_PAIRWISE_KMER_ALGORITHM;
    static string POSITIVE_COMPARISONS_IN_PAIRWISE_KMER_ALGORITHM;
    static string COMPARISONS_IN_SWAT_BRANCH_ALGORITHM;


private:

    class Stat {
    public:
        Stat() {
            maxVal = 1000000000; // 10^9
            maxVal *= -maxVal; // 10^18
            minVal = -maxVal;
            sum = sumOfSquares = numberOfElements = 0;
        }

        double sum;
        double sumOfSquares;
        double numberOfElements;
        double maxVal;
        double minVal;
    };


    static map<string, Stat> data;


};


#endif //GENOMEALIGNMENT_STATISTICSGENERATORBIGDATA_H
