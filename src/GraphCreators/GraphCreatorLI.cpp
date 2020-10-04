//
// Created by sylwester on 12/26/18.
//

#include <GraphCreators/GraphCreatorLI.h>
#include <GraphCreators/GraphCreatorPairwiseKmerBranch.h>

GraphCreatorLI::GraphCreatorLI(vector<Read *> *reads, Graph *G) : GraphCreatorKmerBased(reads, G) {
    graphCreator = new GraphCreatorPairwiseKmerBranch(reads, G);
}


void GraphCreatorLI::createAlignmentsForKmers(vector<Kmer> &kmers, int p, int q, int thread_id) {
    graphCreator->createAlignmentsForKmers(kmers, p, q);
}


void GraphCreatorLI::startAlignmentGraphCreation() {

    for (int i = 0; i < min(4, Params::LI_PRIORITIES_TO_CONSIDER); i++) {
        cerr
                << "******************************************************************************************************************  STARTING ALIGNMENT GRAPH CREATION FOR PRIORITIES:\t";
        WRITE1(Read::priorities);
        graphCreator->startAlignmentGraphCreation();
        rotate(Read::priorities.begin(), Read::priorities.begin() + 1, Read::priorities.end());
    }

}
