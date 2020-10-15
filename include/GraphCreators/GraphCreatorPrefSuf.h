//
// Created by sylwester on 12/31/18.
//

#ifndef GENOMEALIGNMENT_GRAPHCREATORPREFSUF_H
#define GENOMEALIGNMENT_GRAPHCREATORPREFSUF_H


#include "GraphCreator.h"
#include "DataStructures/Kmer.h"
#include<deque>
#include <DataStructures/KmerGCPS.h>
#include <atomic>

/*
 * Cherry
 */
class GraphCreatorPrefSuf : public GraphCreator {
public:
    GraphCreatorPrefSuf(vector<Read *> *reads, Graph *G, bool remove_isolated_reads = false);

    virtual ~GraphCreatorPrefSuf();

    void clear();

    void startAlignmentGraphCreation() override;

    static void test();


private:


    int maxReadLength;
    int currentPrefSufLength;
    Params::KMER_HASH_TYPE prefHashFactor; // this is the hashFactor to create prefix hashes.


    using ADDITIONAL_HASH_TYPE = unsigned;
    const ADDITIONAL_HASH_TYPE MAX_ADDITIONAL_HASH = 1'000'000'007; //  500'009; // is also a prime
    ADDITIONAL_HASH_TYPE prefHashFactorAdditional;

    /*
     * This is an array of prefixes of all reads. prefixKemrs[i] is the hash of current prefix of read with id i
     */
    vector<unsigned long long> prefixKmers; // just the hash values of kmers
    vector<ADDITIONAL_HASH_TYPE> prefixKmersAdditional; // just the additinal hash check

    /**
     * Analogous to prefixKmers
     */
    vector<unsigned long long> suffixKmers;

    vector<ADDITIONAL_HASH_TYPE> suffixKmersAdditional; // just the additional hash check

    /**
     * This is used to keep only Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN edges with small overlap.
     * It will be used only if Params::SPACE_EFFICIENT is set to 1.
     */
    static const int SOES = 3;
//    typedef pair<unsigned, unsigned> *SOES_TYPE;
//    vector<SOES_TYPE> smallOverlapEdges;


    void removeKmersFromBucketsJob(int a, int b, int thread_id);

    void putKmersIntoBucketsJob(int a, int b, int thread_id);

    /*
     * This is the vector that keeps prefixes (ids) in buckets.
     */
    vector<vector<unsigned> > prefixKmersInBuckets; // 'hash map' to store just the ids of reads, from which prefixes are considered
    int prefixKmersBuckets; // this is the number of buckets.

    void calculateMaxReadLength();

    /**
     * Initializes necessary vectors and creates hashes for all suffixes and prefixes of minimum required overlap (perhaps +-1 of that length).
     */
    void createInitialState();


    bool updatePrefixHash(int id, int currentPrefSufLength, Params::KMER_HASH_TYPE prefHashFactor,
                          ADDITIONAL_HASH_TYPE prefHashFactorAdditional); // function updates hash of prefix of read with given id (Moves one position to the right). Returns tru if hash was updated or false if hash is already the length of the read
    bool updateSuffixHash(int id,
                          int currentPrefSufLength); // function updates hash of suffix of read with given id (moves one position to the left). Returns tru if hash was updated or false if hash is already the length of the read




    void nextPrefSufIteration();


    void updatePrexihHashJob(int a, int b, int thread_id);

    void nextPrefSufIterationJobAddEdges(int a, int b, int thread_id);

    void createInitialStateJob(int a, int b, int thread_id);

    void moveSmallOverlapEdgesToGraphJob(int a, int b, int thread_id);

    /**
     * If true, then isolated reads will be removed from graph before reversing the graph. It is to reduce memory peak.
     */
    const bool removeIsolatedReadsBeforeReversingGraph;

    void writeState();


    // STATISTICS
    std::atomic<long long> bitsetChecksCount;
    std::atomic<long long> goodBitsetChecksCount;
    std::atomic<long long> bitsetCheckEdgesRemoved;
    std::atomic<long long> prefSufChecks;
    std::atomic<long long> goodPrefSufChecks;

};


#endif //GENOMEALIGNMENT_GRAPHCREATORPREFSUF_H
