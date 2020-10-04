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

