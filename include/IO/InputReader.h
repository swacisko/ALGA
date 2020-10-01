/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InputReader.h
 * Author: sylwester
 *
 * Created on November 23, 2018, 7:59 PM
 */

#ifndef INPUTREADER_H
#define INPUTREADER_H

#include<string>
#include <future>
#include "StatisticsGenerators/StatisticsGenerator.h"
#include "DataStructures/Read.h"
#include "Global.h"

using namespace std;

class InputReader {
public:
    InputReader();
//    InputReader(const InputReader& orig);
//    virtual ~InputReader();

    void readInput();

    string readOneRead1(istream &str);


private:
    void readReads();

    bool readPairedReads(); // reads all files from file with name *2 (this is a file that contains all paired reads.

    VI STRreads;
    VI Nreads;

    VVI NsInRead;


    void readParallelJob(vector<vector<Read *> > &reads, int thread_id);

};

#endif /* INPUTREADER_H */

