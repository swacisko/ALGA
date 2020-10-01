//
// Created by sylwester on 7/2/19.
//

#ifndef GENOMEALIGNMENT_CONTIGCORRECTOR_H
#define GENOMEALIGNMENT_CONTIGCORRECTOR_H


#include <DataStructures/Graph.h>

class ContigCorrector {
public:

    ContigCorrector(Graph *G, vector<Read *> &reads);


    /**
     * Function performs operations aiming at correcting contigs. I includes a.o. trimming ends to reduce duplication ratio.
     * @param contigs
     */
    void correctContigs(vector<Contig *> &contigs);


    /**
     * Function checks whether contigs have same nodes in the graph. If they have, some of these contigs are trimmed, to eliminate duplication ratio
     */
    void trimContigsWithSameEnds(vector<Contig *> &contigs);


private:

    /**
     *
     * @param contigs
     * @return pair of vectors. first vector contains number of contigs that begin in given nodes, second vector contains number of contigs that end in given node
     */
    pair<VI, VI> createGraphForContigEnds(vector<Contig *> &contigs);

    Graph *G;
    vector<Read *> *reads;


    /**
     *
     * @param ends
     * @return pair of set. First set contains ids of contigs that should have their beginning trimmed, second set contains ids of contigs that should have their end trimmed
     */
    pair<set<int>, set<int> > markContigsToTrim(vector<Contig *> &cotnigs, pair<VI, VI> &ends);

};


#endif //GENOMEALIGNMENT_CONTIGCORRECTOR_H
