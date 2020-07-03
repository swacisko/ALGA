/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Global.cpp
 * Author: sylwester
 * 
 * Created on November 24, 2018, 6:03 PM
 */

#include <Global.h>
#include <Utils/MyUtils.h>
#include <future>

#include "Global.h"

Global::Global() {
}

Global::Global(const Global& orig) {
}

Global::~Global() {
}




void Global::addRead(string s, vector<Read*> & READS) { // adds read to the set of reads.

        Read* r = new Read(READS.size(),s) ;
        READS.push_back( r );
}



vector<Read*> Global::READS(0);
Graph Global::GRAPH(0);
double Global::AVG_READ_LENGTH = 0;

void Global::writeReads(int a, int b) {
    for(int i=a; i<=b; i++){
        cerr << *READS[i] << endl;
    }
}

void Global::removeRead(int id) {
    if( READS[id] != nullptr ){
        delete READS[id];
        READS[id] = nullptr;
    }
}

void Global::removeIsolatedReads() {

    VI* indeg = GRAPH.getInDegrees();
//
//
//    for( int i=0; i<GRAPH.size(); i++ ){
//        if( (*indeg)[i] == 0 && GRAPH[i].size() == 0 ) removeRead(i);
//    }

    vector< std::future<void> > futures(Params::THREADS-1);

    auto worker = [=, &indeg](int a, int b){
        for( int j=a; j<=b; j++ ){
            if( (*indeg)[j] == 0 && GRAPH[j].size() == 0 ) removeRead(j);
        }
    };

    int W = (int) ceil( (double)GRAPH.size() / Params::THREADS);
    for( int i=1; i<Params::THREADS; i++ ){
        int a = i*W;
        int b = min( (i+1)*W-1, (int)GRAPH.size()-1 );
        futures[i-1] = std::async(std::launch::async, worker, a,b );
    }
    worker(0,W-1);
    for( auto & p : futures ) p.get();

    delete indeg;
    indeg = 0;
}

void Global::generateFasta(string filename) {
    cerr << "CAUTION! IT should work if there are REVCOMP reads and PAIRED read, it may not work properly if there are no REVCOMP or PAIRED reads" << endl;

    ofstream fastastream(filename + "_1.fasta", ios::out);

    const string N_STRING = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"; // 100 N's
    int progressCounter = 0;
    for( int i=0; i<READS.size(); i++ ){
        if( Params::ADD_COMP_REV_READS && ( i % 2 == 0 ) ) continue;
        if( Params::ADD_PAIRED_READS && ( i % 4 == 3 ) ) continue;  // this is here to avoid writin piared reads (from file 2)

        fastastream << ">" << i << endl;
        if( READS[i] != nullptr ) fastastream << READS[i]->getSequenceAsString() << endl;
        else fastastream << N_STRING << endl;

        MyUtils::writeProgress( i+1, READS.size(), progressCounter, "creating first file",1 );
    }
    cerr << endl;
    fastastream.close();

    fastastream.open(filename + "_2.fasta", ios::out);

    progressCounter = 0;
    for( int i=0; i<READS.size(); i++ ){
        if( Params::ADD_COMP_REV_READS && ( i % 2 == 0 ) ) continue;
        if( Params::ADD_PAIRED_READS && ( i % 4 == 1 ) ) continue;  // this is here to avoid writin paired reads (from file 1)

        fastastream << ">" << i << endl;
        if( READS[i] != nullptr ) fastastream << READS[i]->getSequenceAsString() << endl;
        else fastastream << N_STRING << endl;

        MyUtils::writeProgress( i+1, READS.size(), progressCounter, "creating second file",1 );
    }

    fastastream.close();
    cerr << endl;

}