//
// Created by sylwester on 12/31/18.
//

#include <GraphCreators/GraphCreator.h>
#include <Global.h>

#include "GraphCreators/GraphCreator.h"

GraphCreator::GraphCreator(vector<Read *> *reads, Graph *G) : reads(reads), G(G) {
    alignFrom = VB(G->size(), true);
    alignTo = VB(G->size(), true);
}

GraphCreator::~GraphCreator() {
    reads = 0;
    G = 0;
}


//void GraphCreator::addCompRevConnection(Read *r1, Read *r2, int offset) {
//    Read * r1revcomp = Global::READS[ r1->getIdOfCompRevRead() ];
//    Read * r2revcomp = Global::READS[ r2->getIdOfCompRevRead() ];
//    int rightOffset = Read::getRightOffset( r1,r2,offset );
//
//    addDirectedEdge( r2revcomp->id, r1revcomp->id, rightOffset );
//}
