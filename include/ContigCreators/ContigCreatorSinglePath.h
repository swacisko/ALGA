//
// Created by sylwester on 3/4/19.
//

#ifndef GENOMEALIGNMENT_CONTIGCREATORSINGLEPATH_H
#define GENOMEALIGNMENT_CONTIGCREATORSINGLEPATH_H


#include <unordered_set>
#include <unordered_map>
#include "ContigCreator.h"

class ContigCreatorSinglePath : public ContigCreator {
public:

    ContigCreatorSinglePath(Graph *G, vector<Read *> &reads);


    vector<Contig *> getAllContigs() override;

    static void test();

private:
    vector<Contig *> getContigOmitShortCyclesFrom(int beg);


    void addContractedPathToString(int a, int b, string &s, vector<pair<Read *, int>> &readsInContig);

    void addContractedPathToString(int a, LPII &path, string &s, vector<pair<Read *, int>> &readsInContig);

    void correctSNPsJob(int a, int b, int thread_id, vector<Contig *> &contigs);

    /**
     * Function returns a vector of possible candidates for the next vertex during graph traverse. It is of the form of pairs (nextId, offset)
     * @param beg
     * @param p
     * @param readsInContig
     * @return
     */
    VPII getNextStepCandidates(int predecessor, int p, vector<pair<Read *, int>> &readsInContig);

    /**
     * Function traverses the graph and all fork (checks all possibilities) up to length maxLengthOfInsertSize from node p. As next step candidates
     * returns those neighbors of p that lie on path containing enough paired connections.
     * @param predecessor
     * @param p
     * @param readsInContig
     * @return
     */
    VPII getNextStepCandidatesByPairedReads(int predecessor, int p, vector<pair<Read *, int>> &readsInContig);

    /**
     *
     * @param p
     * @param d
     * @param readsInContig
     * @return true if node d can be next step candidate of p during traverse
     */
    bool canBeNextStepCandidate(int predecessor, int p, int d, int of,
                                vector<pair<Read *, int>> &readsInContig);





    //*******************//*******************//*******************//*******************//*******************


    /**
     * Function maps all reads that are in edges to those edges.
     */
    void mapReadsToContigs();


    bool canBeNextStepCandidateByPairedReads2(int predecessor, int p, int d, int of,
                                              vector<pair<Read *, int>> &readsInContig);


    //**************************//**************************//**************************//**************************

    /**
     * For each node a with outdeg 1 ( edge (a,b) with length >= 2*avg_read_length), we check all its predecessors, that is all nodes d with (d,a) as an edge. From all suchs d's with length((d,s)) >= 2*avg_read_length
     * we select the one that has greatest number of paired connections in it
     */
    void markReliablePredecessorsByPairedConnections();

    /**
     *
     * @param a
     * @return id of a reliable predecessor of node a, or -1 if a has no reliable predecessors.
     */
    int getReliablePredecessor(int a);

    /**
     * The same as getReliablePredecessor, but we admit multiple predecessors possible.
     * @param a
     * @return
     */
    VI getReliablePredecessors(int a);

    /**
     * Calculated number of pairs of reads (r1,r2) such that r1 and r are paired and r1 is on edge (d,a), r2 is on (a,b)
     * @param d
     * @param a
     * @param b
     * @return
     */
    int countPairedConnections(int d, int a, int b);

    /**
     * reliablePredecessors[i] is the id of the reliable predecessor of i or -1 if i has no such predecessor
     */
//    VI reliablePredecessors;
//    vector< unordered_set<int> > reliablePredecessors;
    unordered_map<int, unordered_set<int> > reliablePredecessors;

    /**
     * Reverse graph of G for calculating reliable predecessors.
     */
//    VVPII GRev;
    unordered_map<unsigned, VPII> GRev;

    /**
     * Only edges of that length will be considered for predecessor check
     */
    int minLengthOfEdgeForReliablePredecessor;

    /**
     * This is not to check all reads on the edge, just the beginning or the end of the edge, of that length, since paired reads cannot be in greater distance
     */
    static const int maxLengthOfInsertSize = 1000;

    static const int minPairedConnections = 5;



    //**************************//**************************//**************************//**************************

    VI *inDeg;
    VB was;
    VI dst;
//    VVPII* initialStructure; // this is the structure of the initial graph. Initial graph is much more dense, than the resulting one,
};


#endif //GENOMEALIGNMENT_CONTIGCREATORSINGLEPATH_H
