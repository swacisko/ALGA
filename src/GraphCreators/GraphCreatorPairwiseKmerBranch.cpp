//
// Created by sylwester on 5/20/19.
//

#include "GraphCreators/GraphCreatorPairwiseKmerBranch.h"


GraphCreatorPairwiseKmerBranch::GraphCreatorPairwiseKmerBranch(vector<Read *> *reads, Graph *G) : GraphCreatorKmerBased(
        reads, G) {

    ach = alignmentController;
    neighbors = VI(G->size(), (int) Params::INF);
}


void GraphCreatorPairwiseKmerBranch::createAlignmentsForKmers(vector<Kmer> &kmers, int p, int q, int thread_id) {

    int D = q - p + 1;

    int LARGE_BUCKET_SIZE = 180;
    if (branchMarkers.size() < D || (branchMarkers.size() > LARGE_BUCKET_SIZE && D < LARGE_BUCKET_SIZE))
        branchMarkers = vector<Bitset>(D, Bitset(D));
    else {
        for (int i = 0; i < D; i++) { // clearing for marking branches
            branchMarkers[i].set(0, D, 0);
        }
    }


    int WINDOW_SEARCH_SIZE = Params::INF;


    for (int i = q - 1; i >= p; i--) {
        int id1 = kmers[i].read->getId();
        if (!alignFrom[id1]) continue;

        int indInRead1 = kmers[i].indInRead;

        G->lockNode(id1);
        for (PII &x : (*G)[id1]) neighbors[x.first] = x.second;
        G->unlockNode(id1);


        for (int j = i + 1; j <= min(q, i + WINDOW_SEARCH_SIZE); j++) {

            int id2 = kmers[j].read->getId();
            if (!alignTo[id2]) continue;
            if (id1 == id2) continue;

            int indInRead2 = kmers[j].indInRead;
            int offset = indInRead1 - indInRead2;
            if (offset < Params::MIN_OFFSET_FOR_ALIGNMENT) continue;


            if (100 * offset > Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT * kmers[i].read->size()) break;


            int overlap = calculateReadOverlap(kmers[i].read, kmers[j].read, offset);
            if (overlap < Params::MIN_OVERLAP_AREA) continue;


            if (Read::getRightOffset(kmers[i].read, kmers[j].read, offset) < 0) continue;


            G->lockNode(id1);

            if (branchMarkers[i - p][j - p] == 0) { // if i cannot yet get from id1 to id2

                if (neighbors[id2] > offset) { // if there is an edge, but with offset greater than current offset

                    if (ach->canAlign(kmers[i].read, kmers[j].read, offset)) {
                        G->addDirectedEdge(id1, id2, offset);
                        neighbors[id2] = offset;
                    }

                }

                if (neighbors[id2] !=
                    Params::INF) { // if there was and edge, but longer than current offset, or i just added it, i mark branchMarkers
                    branchMarkers[i - p].set(j - p, 1);
                    branchMarkers[i - p] |= branchMarkers[j - p];
                }

            }


            G->unlockNode(id1);
        }

        G->lockNode(id1);
        for (PII &x: (*G)[id1]) neighbors[x.first] = Params::INF;
        G->unlockNode(id1);

    }


}

