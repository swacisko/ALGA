//
// Created by sylwester on 5/13/19.
//

#ifndef GENOMEALIGNMENT_GRAPHVISUALIZER_H
#define GENOMEALIGNMENT_GRAPHVISUALIZER_H


#include <DataStructures/Graph.h>

class GraphVisualizer {

public:


    void writeWholeGraph(Graph *G, vector<Read *> &reads, string filename);


    /**
     * Writes the graph to .gv format to use with Graphviz visualization tools.
     * Function creates graph as follows:
     * 1. Only those connected components of this graph are considered that contain at least one contig for @contigs.
     * 2. Only nodes that are a begginning of a contig, and end of a contig or have indegree or outdegree 0 or >= 2 will be drawn as nodes, other will be in contracted paths
     * 3. Edges will have weight equal to the length of the contracted path they represent.
     *
     */
    void writeInGraphvizFormat(Graph *G, vector<Contig *> &contigs);


private:

    set<int> relevantNodes;
    ofstream out;

    static const int C = 33;
    string colors[C] = {"blue", "maroon", "green", "orange", "red", "brown", "beige", "cyan", "gold", "darkgreen",
                        "darkorchid", "goldenrod4", "gray",
                        "aquamarine", "pink", "deeppink", "coral", "cornflowerblue", "indigo", "sienna", "turquoise",
                        "salmon", "coral3", "forestgreen", "greenyellow",
                        "hotpink", "magenta1", "mediumblue", "plum", /*"snow",*/ "yellow3", "wheat", "turquoise",
                        "crimson"};


    void createRelevantNodes(Graph *G, vector<Contig *> &contigs);

    void visualizeGraph(Graph *G, vector<Contig *> &contigs);

    void visualizeContigs(vector<Contig *> &contigs);

    void visualizeContig(Contig *ctg, string color);


    bool ADD_REVCOMP_CONTIGS;
};


#endif //GENOMEALIGNMENT_GRAPHVISUALIZER_H
