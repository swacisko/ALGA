//
// Created by sylwester on 5/20/19.
//

#ifndef GENOMEALIGNMENT_GRAPHCREATORPAIRWISEKMERBRANCH_H
#define GENOMEALIGNMENT_GRAPHCREATORPAIRWISEKMERBRANCH_H


#include "GraphCreatorKmerBased.h"

/**
 * This is version of PairwiseKmer GraphCreator that does not create full graph but just branch paths like SwatBranch. This creator however uses exhaustive check and marks on bitsets nodes to which one can move.
 */
class GraphCreatorPairwiseKmerBranch : public GraphCreatorKmerBased {
public:
    GraphCreatorPairwiseKmerBranch(vector<Read *> *read, Graph *G);

    GraphCreatorPairwiseKmerBranch *clone() override { return new GraphCreatorPairwiseKmerBranch(reads, G); };

protected:
    void createAlignmentsForKmers(vector<Kmer> &kmers, int p, int q, int thread_id = 0) override;

private:

    AlignmentController *ach;

    vector<Bitset> branchMarkers;
    VI neighbors;

};


#endif //GENOMEALIGNMENT_GRAPHCREATORPAIRWISEKMERBRANCH_H
