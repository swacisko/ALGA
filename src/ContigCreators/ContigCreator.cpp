//
// Created by sylwester on 3/4/19.
//

#include <ContigCreators/ContigCreator.h>

#include "ContigCreators/ContigCreator.h"

ContigCreator::ContigCreator(Graph *G, vector<Read *> *reads) {
    this->G = G;
    this->reads = reads;


}

ContigCreator::~ContigCreator() {
    reads = 0;
    G = 0;
}
