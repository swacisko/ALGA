/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Graph.h
 * Author: sylwester
 *
 * Created on November 17, 2018, 3:58 PM
 */

#ifndef GRAPH_H
#define GRAPH_H

#include<iostream>
#include<algorithm>
#include<iomanip>
using namespace std;
#include<vector>
#include<map>
#include "Params.h"
#include <list>
//#include<thread>
#include <mutex>
#include <set>
#include "Read.h"
#include "Contig.h"

//#define REP(i,N) for(int i=0; i<N; i++)
typedef vector<int> VI;
typedef vector<VI> VVI;
typedef vector<short> VS;
typedef vector<VS> VVS;

typedef pair<int,int> PII;
typedef vector<PII> VPII;
typedef vector<VPII> VVPII;

typedef list<PII> LPII;
typedef map<int, LPII> MILPII;
typedef vector<MILPII*> VMILPII;



class Graph {
public:
    Graph(int N);
    Graph(const Graph& orig);
    virtual ~Graph();
    
    VPII & operator[](int a){ return V[a]; }
    vector<std::vector<std::pair<int, int>, std::allocator<std::pair<int, int>>>, std::allocator<std::vector<std::pair<int, int>, std::allocator<std::pair<int, int>>>>>::iterator
    begin(){ return V.begin(); }
    vector<std::vector<std::pair<int, int>, std::allocator<std::pair<int, int>>>, std::allocator<std::vector<std::pair<int, int>, std::allocator<std::pair<int, int>>>>>::iterator
    end(){ return V.end(); }
    int size(){ return (int) V.size(); }
    void clear();
    void push_back( /*VPII row*/ ){ /*V.push_back(row);*/ push_node(); }
    void push_node();
//    void pop_back(){ /*V.pop_back();*/ V.pop_back(); /*maps.pop_back();*/ contractedEdges.pop_back(); mutexes->pop_back(); }
    void addDirectedEdge( int a, int b, int offset ); // adds directed edge from a to b with given offset. If such edge exists it will modify the offset to the smallest of thos considered yet for those edge.

    /**
     *
     * @return structure of the graph, vector V.
     */
    VVPII getV(){ return V; }

    /**
     * CAUTION!! This function pushes directed edge (a,b) with weight c to the graph. It does not check whether another edge (a,b) already exists there.
     * It may cause the graph become a multigraph. If we use this function and want to retain only edge with smallest offset, after adding all edges with this method
     * retainOnlySmallestOffset method should be invoked.
     * The graph will still keep only one edge if the node is stored on map.
     *
     *
     * If there are multiple edges
     * @param a
     * @param b
     * @param offset
     */
    void pushDirectedEdge(int a, int b, int offset);

    void removeDirectedEdge( int a, int b ); // removes edge (a,b) from the graph if exists. Changes sequence of neighbors of a during execution.
    void mergeVertices( int a, int b, int offset ); // merges vertices a and b. All edges (a,x) will be added as (b,x). (b,a) will be removed if exists to avoid cycles. (a,b) will be present with given offset.

    /**
     * Function removes all edges that are in given vector.
     * CAUTION! Vector is passed by reference and is modified! (sorted to quickly remove elements)
     * CAUTION! Changes order of elements in @V (sorted to quickly remove
     * @param edges
     */
    void removeDirectedEdges( VPII & edges );

    bool containsEdge( int a, int b );
    bool containsEdgeShorterOrEqual( int a, int b, int offset );
    bool containsEdgeLongerOrEqual(int a, int b, int offset);

    void write();
    /**
     * Writes all nodes that have indeg > 0 or outdeg > 0 with a <= id <= b
     */
    void writeNonisolatedNodes(int a, int b);

    /**
     * Writes all nonisolated nodes
     */
    void writeNonisolatedNodes(){ writeNonisolatedNodes(0, size()-1); }
    void write( int a, int b ); // writes all nodes from given interval (used to debug)

    /**
     *
     * @param i
     * @param k
     * @return V[i][k].second
     */
    int getWeight(int i, int k);

    /**
     *
     * @param a
     * @param b
     * @return weight of edge a->b if exists or -1 otherwise
     */
    int findWeight( int a, int b );

    void clearNode( int v ); // removes all edges that begin in given node
    void clearNeighborsForNode( int v ); // removes all neighbors of given node
    
    VI* getInDegrees();

    /**
     *
     * @return reverse graph to the given one. It contains only basic connections (edges in V or maps vectors, contracted edge will not be taken into account).
     */
    Graph getReverseGraph();

    /**
     * Function reverses the graph. Is result is the same as *this = this->getReverseGraph(), but uses less memory
     *
     * CAUTION, FOR NOW IT IS NOT IMPLEMENTED PROPERLY YET!!
     */
    void reverseGraph();

    /**
     * Performs operations of type vector<PII>().swap( V[i] ) for all nodes
     */
    void pruneGraph();


    vector<pair<int, int>> getNeighbors(int id); // returns all neighbors of given node in form of pairs (neighbor_id, offset)
    
    friend ostream& operator<<( ostream & str, Graph & G );



    void writeAllEdgesWithNode( int id ); // writes all edges that begin or end in given node

    int getGraphNonmapSizeThreshold() const;
//    void setGraphNonmapSizeThreshold(int graphNonmapSizeThreshold);

//    void transformMaps();
//    void mapNodes();

    void serializeGraph(string fileName);
    bool deserializeGraph(string fileName);

    void addCompRevConnection( Read *r1, Read *r2, int offset ); // if r1 align with r2, then r2_revcomp aligns with r1_revcomp with rightOffset. This function adds that connection

    bool operator==( Graph& oth );
    /**
     *
     * @param oth
     * @return true if this graph is a subgraph of oth (i.e. every edge from this graph is contained in the other graph).
     */
    bool operator<( Graph & oth );

    /**
     *  if there are edges (a,b):c and (a,b):d, then only one edge will be retained (a,b):e, where e is minimal weight, e = min(c,d)
     *  CAUTION! All edges must be in V (not in maps) when this function is invoked. Edges in maps will not be considered.
     */
    void retainOnlySmallestOffset(); //



    /**
     * Function sorts all edges in V for all nodes such that V[0].second < V[1].second < ...
     */
    void sortEdgesByIncreasingOffset();



    /**
     * If there is a path a -> b -> c (where indeg(b) = outdeg(b) = 1), then i will contract that path to a -> c and store b in contractedEdges list.
     * @param a
     * @param b
     */
    bool contractPath(int a, int b, int c);

    /**
     *
     * @param a
     * @param b
     * @return  length_of_path / number_of_reads_in_path. This is avg offset on that path
     */
    double getContractedPathDensity(int a, int b);

    void createContractedEdgesVector();

    /**
     *
     * @param a
     * @param b
     * @return a list representing the contracting path. It is of the form (neighbor,offset), so the first node on the path (node a) is not in that list.
     */
    LPII& getContractedEdgePath(int a, int b);
//    LPII empty,temp;

    /**
     * Function writes contracted path from a to b
     * @param a
     * @param b
     */
    void writeContractedPath( int a, int b );

    void lockNode(int id){
        /*cerr << "locking " << id << " / " << mutexes->size() << endl;*/
        (*mutexes)[id].lock();
        /*cerr << "locked" << endl; */
    }
    void unlockNode(int id){
        /*cerr << "unlocking " << id << endl;*/
        (*mutexes)[id].unlock();
        /*cerr << "unlocked" << endl;*/
    }
    vector<mutex> *mutexes;

    void operator=(const Graph &oth);

    /**
     * For each node function removes @perc percent of edges with largest offset
     * @param perc
     */
    void removeEdgesWithLargestOffset( int perc );

    /**
     *
     * @return number of edges c=in the graph
     */
    LL countEdges();

    /**
     * For each edge (a,b):c with weight c < 0 we delete that edge and add edge (b,a):-c
     * @return number of reversed edges.
     */
    int reverseEdgesWithNegativeOffset();

private:


    LL countEdgesJob( int a, int b, int thread_id );

    /**
     * Helper function for parallel computation
     * @param a
     * @param b
     * @param thread_id
     */
    void retainOnlySmallestOffsetJob(int a, int b, int thread_id);


    /**
     * Helper function for parallel computation
     * @param a
     * @param b
     * @param thread_id
     */
    void sortEdgesByIncreasingOffsetJob(int a, int b, int thread_id);


    /**
     * Helper function for parallel creation of reverse graph.
     * @param a
     * @param b
     * @param thread_id
     * @param mutexes
     */
    void getReverseGraphJob(int a, int b, int thread_id, Graph &GRev);


    /**
     * contractedEdges[i] is a map that contains id of neighbo as key and list of nodes on the contracted edge from i to that neighbor.
     */
    VMILPII contractedEdges;

    vector<LPII*> contractedEdgeDummy;

    VVPII V;

    int edges;
};

#endif /* GRAPH_H */
