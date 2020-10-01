/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   GenomeStatisticsCollector.cpp
 * Author: sylwester
 *
 * Created on November 27, 2018, 11:43 PM
 */

#include <valarray>
#include <StatisticsGenerators/GenomeStatisticsCollector.h>
#include <Global.h>
#include <numeric>
#include <AlignmentControllers/AlignmentControllerLowErrorRate.h>
#include <random>
#include <AlignmentControllers/AlignmentControllerLCS.h>
#include "DataStructures/FAU.h"

#include "StatisticsGenerators/GenomeStatisticsCollector.h"


GenomeStatisticsCollector::GenomeStatisticsCollector() {
}

GenomeStatisticsCollector::GenomeStatisticsCollector(const GenomeStatisticsCollector &orig) {
}

GenomeStatisticsCollector::~GenomeStatisticsCollector() {
}


void GenomeStatisticsCollector::addData(string option, double value) {
//    if( Params::WRITE_STATISTICS == false ) return;
    stats[option] = value;
}

void GenomeStatisticsCollector::writeTestStatistics() {
    cerr << endl << "TEST SPECIFIC STATISTICS:" << endl;
    for (auto a : stats) {
        cerr << a.first << " -> " << a.second << endl;
    }
}


map<string, double> GenomeStatisticsCollector::stats;