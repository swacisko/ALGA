//
// Created by sylwester on 3/4/19.
//

#ifndef GENOMEALIGNMENT_CONTIGCREATOR_H
#define GENOMEALIGNMENT_CONTIGCREATOR_H


#include <DataStructures/Graph.h>
#include <DataStructures/Contig.h>

class ContigCreator {

public:

    ContigCreator(Graph *G, vector<Read *> *reads);

    virtual ~ContigCreator();

    virtual vector<Contig *> getAllContigs() = 0;


    VLL getContigsLengths() { return contigsLengths; }

protected:
    Graph *G;
    vector<Read *> *reads;

    VLL contigsLengths;
};


#endif //GENOMEALIGNMENT_CONTIGCREATOR_H
