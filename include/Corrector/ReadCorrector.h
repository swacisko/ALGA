//
// Created by sylwester on 2/6/20.
//

#ifndef ALGA_READCORRECTOR_H
#define ALGA_READCORRECTOR_H

#include "DataStructures/Read.h"

#include <unordered_map>


class ReadCorrector {
public:
    typedef LL BIG_TYPE;
    typedef int SMALL_TYPE;
    typedef unordered_map <BIG_TYPE, unordered_map<SMALL_TYPE, int>> MAP_TYPE;

    ReadCorrector(vector<Read *> &reads, int sLength, int bLength);

    void correct();

    static void test();

private:

    vector<Read *> *reads;

    int smallLength, bigLength;

    /**
     * If true, then all reads will be reversed and corrected form other side
     */
    bool reversedReads;

    static const int candidateThreshold = 2;
    static const int FREQS_SIZE = (1 << 20);
    vector<mutex> mutexes;


    /**
     * frequencies[x] is a map that contains numbers of occurences of short sequences before sequence x.
     * That is each kmer of length
     */
    vector<MAP_TYPE> frequencies;


    void correctInDirection(bool reversed = false);

    /**
     * Creates @{frequencies}
     */
    void createFrequenciesMap();

    void addReadDataToMap(Read *r, vector<MAP_TYPE> &map);


    void applyCorrection();

    void applyCorrectionToRead(Read *r);


    /**
     * Helper function to avoid actually reversing reads, but just indexing them backwards.
     * @param r
     * @param pos
     * @return
     */
    int accessReadPosition(Read *r, int pos);

    void setReadAtPosition(Read *r, int pos, int val);

    void debugFrequencies();

    string hashToString(LL val, int l);

};

#endif //ALGA_READCORRECTOR_H
