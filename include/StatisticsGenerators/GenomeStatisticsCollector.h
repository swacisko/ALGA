/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GenomeStatisticsCollector.h
 * Author: sylwester
 *
 * Created on November 27, 2018, 11:43 PM
 */

#ifndef GENOMESTATISTICSCOLLECTOR_H
#define GENOMESTATISTICSCOLLECTOR_H

#include<vector>
#include "Utils/MyUtils.h"
#include "IO/InputReader.h"


#include<map>
#include <random>

using namespace std;
typedef long long LL;

class GenomeStatisticsCollector {
public:
    GenomeStatisticsCollector();

    GenomeStatisticsCollector(const GenomeStatisticsCollector &orig);

    virtual ~GenomeStatisticsCollector();

    static void addData(string option, double value);

    static void writeTestStatistics();


private:
    static map<string, double> stats;
};

#endif /* GENOMESTATISTICSCOLLECTOR_H */

