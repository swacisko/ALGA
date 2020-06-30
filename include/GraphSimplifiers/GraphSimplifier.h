/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GraphSimplifier.h
 * Author: sylwester
 *
 * Created on November 25, 2018, 12:42 AM
 */

#ifndef GRAPHSIMPLIFIER_H
#define GRAPHSIMPLIFIER_H

#include "DataStructures/Graph.h"
#include<map>
#include<set>
#include <DataStructures/FAU.h>

typedef vector<VB> VVB;

class GraphSimplifier {
public:
    GraphSimplifier(Graph &G, vector<Read *> &reads);
    GraphSimplifier(const GraphSimplifier& orig);
    virtual ~GraphSimplifier();
    void clear();
    
    /*
     Function removes from graph G every edge (a,c) such that there exist edges (a,b) and (b,c) with property W[a][b] + W[b][c] <= W[a][c]
     */
    void cutNonAndWeaklyMetricTriangles();



    /**
     * If there is en edge with offset 0, i merge the shorter read into the longer one.
     */
    int mergeLength0Edges();

    void simplifyGraphOld(); // function performs several operations (algorithms) to prune read graph of unneccessary edges.
    void simplifyGraph(); // function performs several operations (algorithms) to prune read graph of unneccessary edges.

    /*
        removes short paths that run paralelly. e.g a->b->c and a->d->c are parallel, but will not be removed in removeMetricTriangles. This may be used instead of cutNonAndWeaklyMetricTriangles
     * because id does much more (and it also removes metric triangles),
     */
    void removeShortParallelPaths(int maxOffset); //


    /**
     * Function contracts all paths in graph to one long path, as long as it is possible.
     * @return true if there was at least one contracion done. False if no contraction was done
     */
    bool contractPathNodes();



    /**
     * Function retains in graph only those edges that represent aligning pair of reads with overlap area at least @min_overlap_to_retain and @number_of_long_edges_to_retain edges among other edges,
     * that have the greatest values of overlap area. E.g if there are edges (a,b):50(ov), (a,c):40(ov) and (a,d):30(ov) and @min_overlap_to_retain == 45 and @number_of_long_edges_to_retain == 1 then only edges (a,b):50, (a,c):40
     * will be retained in the graph (where notation (a,b):c(ov) denotes directed edge a->b with corresponding overlap c)
     * @param min_overlap_to_retain
     * @param number_of_long_edges_to_retain
     */
    void removeSmallOverlapEdges(int min_overlap_to_retain, int number_of_long_edges_to_retain);


    /**
     * For each short edge (a,b) we check N(b) in graph G and N(a) in graph GRev.
     * Ich in those neighborhoods there are no paired rads, then we remove that edge.
     */
    bool removeEdgePairedEndSeparators();


    static void test();

private:
    Graph * G;
    vector<Read*> * reads;
//    VI *inDeg;

    /**
     * edgesToRemove[i] is the vector of edges that are to be removed that were added by i-th thread. So edgesToRemove.size() == Params::THREADS
     */
    VVPII edgesToRemove;

//    VVI dst;
//    VVI par;
    VVB was;

    

    void tryToRemoveShortPathsMST(int beg, int maxOffset, int thread_id);


    /*
     * Function removes from graph G dangling branches.
     * A dangling branch is a path starting from node with outdeg >= 2 and ending in node with outdeg = 0. This path contains no vertices (other than the beginning) with outdeg >= 2.
     * Such a path will be removed if it has total length <= maxOffset.
     */
    int removeDanglingBranches(int maxOffset);



    /*
     * It is helper for removeDanglingBranches. Removes branches for given beginning node (adds edges that should be removed to edgesToRemove vector, they will be removed simulataneously in
     * removeDanglingBranches function).
     */
    int removeDanglingBranchesFromNode(int beg, int maxOffset, vector<pair<int, int> > &edgesToRemove,
                                       int thread_id);

    /*
     * Removes dangling branches that start in node with indegree = 0 and end in node with indegree >= 2. This is to ensure that we have small repetitions in regions covered by contigs.
     * This however will not reduce the number of contigs (because removed dangling upper branches will not reduce outdegree of vertices with outdeg >= 2 and will not be produces as contigs
     * for their length will be less than Params::CONTIG_MIN_OUTPUT_LENGTH ).
     */
    int removeDanglingUpperBranches(int maxOffset);


    /**
     * Helper function for parallel execution in cutting metric triangles.
     * @param a
     * @param b
     * @param GRev
     */
    void removeNonAndWeaklyMetricTrianglesJobAddToRemove(int a, int b, int thread_id);

    /**
     * Helper function used to parallelly remove edges that are to be removed in removeNonAndWeaklyMetricTriangles.
     * @param a
     * @param b
     * @param thread_id
     */
    void removeNonAndWeaklyMetricTrianglesJobRemoveEdges( int a, int b, int thread_id );

    void removeShortParallelPathsJob(int a, int b, int maxOffset, int thread_id);

    /**
     * Merges node a into node b
     * @param a
     * @param b
     */
    void mergeNodes(int a, int b, Graph &GRev, int offsetAB);


    void removeSmallOverlapEdgesJob(int a, int b, int min_overlap_to_retain, int number_of_long_edges_to_retain,
                                    int thread_id);


    void clearEdgesToRemove();

    void removeEdgesToRemove();

};

#endif /* GRAPHTRIANGLEREMOVER_H*/

