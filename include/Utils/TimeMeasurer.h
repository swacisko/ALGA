/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TimeMeasurer.h
 * Author: sylwester
 *
 * Created on November 28, 2018, 4:25 PM
 */

#ifndef TIMEMEASURER_H
#define TIMEMEASURER_H

#include<iostream>
#include<map>

using namespace std;
typedef long long LL;
typedef pair<LL,LL> PLL;

#include<ctime>

class TimeMeasurer {
public:
    TimeMeasurer();
    TimeMeasurer(const TimeMeasurer& orig);
    virtual ~TimeMeasurer();
    
    static void stopMeasurement( string option );
    static void startMeasurement( string option );
    static float getMeasurementTimeInSeconds(string option);
    static map<string,float> getAllMeasurements(); // returns all measurements in second
    static void writeAllMeasurements(); // writes all measurements
    static void clearOption( string option ); // clears space in times for given option.
    
    static string INPUT_READER;
    static string GRAPH_SIMPLIFIER;
    static string GRAPH_CREATOR;
    static string ALIGNMENT_CONTROLLER_CAN_ALIGN_BITMAP;
    static string ALIGNMENT_CONTROLLER_CAN_ALIGN_LCS;
    static string OUTPUT_WRITER;
    static string TOTAL_TIME;
    static string GRAPH_LCS_CHECK_BEFORE_REMOVING_PARALLEL_PATHS;
    static string GRAPH_LCS_CHECK_AFTER_REMOVING_PARALLEL_PATHS;
    static string KMER_BUCKETS_SORTING;

    
private:

    
    
    
    static map<string,LL> times;
    static map<string,LL> timesTotal; // total time of measurement for given parameter in clock() units (CLICKS_PER_SEC). Sum of times between startMEasurement() and stopMeasurement()
};

#endif /* TIMEMEASURER_H */

