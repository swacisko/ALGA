//
// Created by sylwester on 12/26/18.
//

#ifndef GENOMEALIGNMENT_GRAPHCREATORLI_H
#define GENOMEALIGNMENT_GRAPHCREATORLI_H


#include "GraphCreatorKmerBased.h"

class GraphCreatorLI : public GraphCreatorKmerBased {

public:
    GraphCreatorLI(vector<Read *> *reads, Graph *G);

    void startAlignmentGraphCreation() override;

    void createAlignmentsForKmers(vector<Kmer> &kmers, int p, int q, int thread_id = 0) override;

    GraphCreatorLI *clone() override { return new GraphCreatorLI(reads, G); }

    void setAlignFrom(int from, bool val) override {
        alignFrom[from] = val;
        graphCreator->setAlignFrom(from, val);
    }

    void setAlignTo(int to, bool val) override {
        alignTo[to] = val;
        graphCreator->setAlignTo(to, val);
    }

protected:

//    vector< vector<Kmer> > getKmersForBucket(int bucket) override;

    GraphCreatorKmerBased *graphCreator;
};


#endif //GENOMEALIGNMENT_GRAPHCREATORLI_H
