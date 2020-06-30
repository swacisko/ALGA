/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   StatisticGenerator.h
 * Author: sylwester
 *
 * Created on November 27, 2018, 1:30 PM
 */

#ifndef STATISTICGENERATOR_H
#define STATISTICGENERATOR_H

#include<iostream>
#include<vector>
#include<cmath>
#include<map>
#include "Params.h"
//#include<algorithm>
#include<numeric>

using namespace std;

/*
 CAUTION. All methods may change the order of elements in given vector of data.
 */
class StatisticsGenerator {
public:
    StatisticsGenerator();
    StatisticsGenerator(const StatisticsGenerator& orig);
    virtual ~StatisticsGenerator();
    
    template<class _T>
    static double getAverage(vector<_T> data) {
        sort( data.begin(), data.end() );
        double s = getSum(data);
        if( data.size() == 0 ) return 0;
        return s / data.size();
    }

    template<class _T>
    static _T getMax( vector<_T> data ) {
        return *max_element(data.begin(), data.end());
    }

    template<class _T>
    static _T getMedian(vector<_T> data) {
        sort( data.begin(), data.end() );
        return data[ data.size() >> 1 ];
    }

    template<class _T>
    static _T getMin(vector<_T> data) {
        return *min_element(data.begin(), data.end());
    }

    template<class _T>
    static double getStandardDeviation(vector<_T> data) {
        double avg = getAverage(data);
        double stddev = 0;
        for( auto a : data ){
            double v = a-avg;
            v *= v;
            stddev += v;
        }

        stddev /= data.size();
        stddev = sqrt(stddev);
        return stddev;
    }

    template<class _T>
    static map<string, double> generateAllStatistics(vector<_T> data) {
        sort( data.begin(), data.end() );
        map<string, double> stats;
        if( data.empty() ) return stats;

        stats["number of elements"] = (double)data.size();
        stats[ "sum of elements" ] = (double) getSum(data);
        stats[  "average" ] = (double) getAverage( data  );
        stats[ "minVal element"  ] = (double) getMin( data );
        stats[ "maxVal element"  ] = (double) getMax( data );
        stats[ "median"  ] = (double) getMedian( data );
        stats[ "standard deviation"  ] = (double) getStandardDeviation(data);
        stats[  "sum of squares of elements" ] = (double) getSumOfSquares( data );
        return stats;
    }

    template<class _T>
    static void writeAllStatistics(vector<_T> data) {
//        if( Params::WRITE_STATISTICS == false ) return;

        sort( data.begin(), data.end() );
        map<string,double>  stats = generateAllStatistics(data);
        for(auto a : stats ) cerr << a.first << "  -->  " << a.second << endl;
        cerr << endl;
    }

    template<class _T>
    static void writeAllStatistics( map<_T,_T> data ) {
//        if( Params::WRITE_STATISTICS == false ) return;

        map<string,double> stats;
        stats["number of elements"] = accumulate( data.begin(), data.end(), 0ll, [](auto a, auto b){ return a + b.second; } );
        stats[ "sum of elements" ] =  accumulate( data.begin(), data.end(), 0ll, [](auto a, auto b){ return a + b.first * b.second; } );
        stats[ "average" ] = stats[ "sum of elements" ] / stats[ "number of elements" ];
        stats[ "minVal element"  ] = accumulate( data.begin(), data.end(), 1000000000ll, [](auto a, auto b){ return min( a, b.first ); } );
        stats[ "maxVal element"  ] = accumulate( data.begin(), data.end(), 0ll, [](auto a, auto b){ return max( a, b.first ); } );

        double avg = stats["average"];
        stats[ "standard deviation"  ] = sqrt( (double)accumulate( data.begin(), data.end(), 0ll, [&avg](auto a, auto b){
                                                                    return a + b.second * ( b.first-avg ) * (b.first - avg); } ) / stats["number of elements"] );
        stats[  "sum of squares of elements" ] = accumulate( data.begin(), data.end(), 0ll, [](auto a, auto b){ return a + b.second * b.first * b.first; } );

//        stats["number of elements"] = (double)data.size();
//        stats[ "sum of elements" ] = (double) getSum(data);
//        stats[  "average" ] = (double) getAverage( data  );
//        stats[ "minVal element"  ] = (double) getMin( data );
//        stats[ "maxVal element"  ] = (double) getMax( data );
//        stats[ "median"  ] = (double) getMedian( data );
//        stats[ "standard deviation"  ] = (double) getStandardDeviation(data);
//        stats[  "sum of squares of elements" ] = (double) getSumOfSquares( data );

        for(auto a : stats ) cerr << a.first << "  -->  " << a.second << endl;
    }

    template<class _T>
    static vector<_T> trimEnds(vector<_T> data, int left, int right) {
        vector<_T> v( data.begin()+left, data.begin() + data.size() - right );
        return v;
    }

    template<class _T>
    static _T getSum(vector<_T>& data) {
        _T s = 0;
        for(auto a : data) s += a;
        return s;
    }

    template<class _T>
    static vector<pair<_T, int> > getCounts(vector<_T>& data) {
        map<_T,int> c;
        for(auto a : data) c[a]++;
        vector< pair<_T,int> > res;
        for(auto a : c){
            res.push_back( {a.first, a.second} );
        }
        return res;
    }

    template<class _T>
    static _T getSumOfSquares(vector<_T>& data) {
        _T s = 0;
        for(auto a : data) s += (double)a*(double)a;
        return s;
    }



    template<class _T>
    static vector<_T> trimEndsInPercent(vector<_T> data, float left, float right) {
        int l = left * data.size();
        int r = right*data.size();
        return trimEnds( data, l,r );
    }

    
private:

};

#endif /* STATISTICGENERATOR_H */

