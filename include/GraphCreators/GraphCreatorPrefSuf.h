//
// Created by sylwester on 12/31/18.
//

#ifndef GENOMEALIGNMENT_GRAPHCREATORPREFSUF_H
#define GENOMEALIGNMENT_GRAPHCREATORPREFSUF_H


#include "GraphCreator.h"
#include "DataStructures/Kmer.h"
#include<deque>

/*
 * Cherry
 */
class GraphCreatorPrefSuf : public GraphCreator {
public:
    GraphCreatorPrefSuf( vector<Read*> *reads, Graph *G );

    virtual ~GraphCreatorPrefSuf();

    void clear();

    void startAlignmentGraphCreation() override;

    static void test();

    void removeSmallOverlapEdges();

private:



    void writeState();

    int maxReadLength;

    /*
     * This is an array of prefixes of all reads. prefixKemrs[i] is the current prefix of read with id i
     */
    vector<Kmer> prefixKmers;

    vector<Kmer> suffixKmers;

    /**
     * This is used to keep only Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN edges with small overlap.
     * It will be used only if Params::SPACE_EFFICIENT is set to 1.
     */
//    VVPII smallOverlapEdges;
    static const int SOES = 3;
    vector< pair<unsigned,unsigned>[SOES] > smallOverlapEdges;

    /**
     * outdegOverThreshold[i] is the number of edges that would begin in node i in the full-errorless-graph, that would have overlap greater than the threshold of 50 * (1 + scale).
     * This is to remove situations where there is a path a->b->c->...-> v and a->w with offset(a->w) below the threshold.
     */
    vector<char> outdegOverThreshold;

    /*
     * prefixModulatorHash[i] is the hash % BUCKET_MODULATOR value.
     */
//    vector< long long > prefixModulatorHash;

    void removeKmersFromBucketsJob(int a, int b, int thread_id);

    void putKmersIntoBucketsJob(int a, int b, int thread_id);

    /*
     * suffixModulator[i] is the hash of i-th suffix modulo BUCKET_MODULATOR value.
     */
//    vector< long long > suffixModulatorHash;
    /*
     * This is the vector that keeps prefixes in buckets. prefixes from prefixKmers are kept in these buckets as pointers.
     */
    vector< vector< Kmer* > > prefixKmersInBuckets;
    int prefixKmersBuckets; // this is the number of buckets.

    void calculateMaxReadLength();

    void createInitialState();
    int currentPrefSufLength;
    Params::KMER_HASH_TYPE prefHashFactor; // this is the hashFactor to create prefix hashes.


    bool updatePrefixHash(int id, int currentPrefSufLength, Params::KMER_HASH_TYPE prefHashFactor); // function updates hash of prefix of read with given id (Moves one position to the right). Returns tru if hash was updated or false if hash is already the length of the read
    bool updateSuffixHash(int id, int currentPrefSufLength); // function updates hash of suffix of read with given id (moves one position to the left). Returns tru if hash was updated or false if hash is already the length of the read




    void nextPrefSufIteration();


    void updatePrexihHashJob( int a, int b, int thread_id );

    void nextPrefSufIterationJobAddEdges( int a, int b, int thread_id );

    void createInitialStateJob( int a, int b, int thread_id );

    void moveSmallOverlapEdgesToGraphJob(int a, int b, int thread_id);

};


#endif //GENOMEALIGNMENT_GRAPHCREATORPREFSUF_H
