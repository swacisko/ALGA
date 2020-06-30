/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TimeMeasurer.cpp
 * Author: sylwester
 * 
 * Created on November 28, 2018, 4:25 PM
 */

#include <Params.h>
#include "Utils/TimeMeasurer.h"

TimeMeasurer::TimeMeasurer() {
}

TimeMeasurer::TimeMeasurer(const TimeMeasurer& orig) {
}

TimeMeasurer::~TimeMeasurer() {
}

void TimeMeasurer::stopMeasurement(string option) {
//    if( Params::WRITE_STATISTICS == false ) return;

    if( times.find(option) == times.end() ) times[option] = clock();
    else{
        timesTotal[option] += ( clock() - times[option] );
        clearOption(option);
    }
}


void TimeMeasurer::startMeasurement(string option) {
//    if( Params::WRITE_STATISTICS == false ) return;
    times[option] = clock();
}

float TimeMeasurer::getMeasurementTimeInSeconds(string option) {
    if( timesTotal.find(option) == timesTotal.end() ) return -1;
    else return ( ( (double)timesTotal[option] / (double)CLOCKS_PER_SEC ) );
}

map<string, float> TimeMeasurer::getAllMeasurements() {
    map<string,float> res;
    for( auto a : timesTotal ){
        res[a.first] = getMeasurementTimeInSeconds(a.first);
    }
    return res;
}

void TimeMeasurer::writeAllMeasurements() {
//    if( Params::WRITE_STATISTICS == false ) return;

    cerr << endl << "TIME MEASUREMENTS:" << endl;
    auto a = getAllMeasurements();
    cerr << fixed;
    cerr.precision(3);
    for( auto b : a ){
        cerr << b.first << " -> " << b.second << " seconds" << endl;
    }
}

void TimeMeasurer::clearOption(string option) {
    times.erase(option);
}

map<string,LL> TimeMeasurer::times;
map<string,LL> TimeMeasurer::timesTotal;
string TimeMeasurer::INPUT_READER = "InputReader";
string TimeMeasurer::GRAPH_SIMPLIFIER = "GraphSimplifier";
string TimeMeasurer::GRAPH_CREATOR = "GraphCreator";
string TimeMeasurer::ALIGNMENT_CONTROLLER_CAN_ALIGN_BITMAP = "AlignmentController canAlign Bitmap";
string TimeMeasurer::ALIGNMENT_CONTROLLER_CAN_ALIGN_LCS = "AlignmentController canAlign LCS";
string TimeMeasurer::OUTPUT_WRITER = "OutputWriter";
string TimeMeasurer::TOTAL_TIME = "TotalTime";
string TimeMeasurer::GRAPH_LCS_CHECK_BEFORE_REMOVING_PARALLEL_PATHS = "Graph LCS check all edges before removing parallel paths";
string TimeMeasurer::GRAPH_LCS_CHECK_AFTER_REMOVING_PARALLEL_PATHS = "Graph LCS check all edges after removing parallel paths";
string TimeMeasurer::KMER_BUCKETS_SORTING = "K-mer bucket sorting";