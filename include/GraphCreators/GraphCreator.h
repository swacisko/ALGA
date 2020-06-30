//
// Created by sylwester on 12/31/18.
//

#ifndef GENOMEALIGNMENT_GRAPHCREATOR_H
#define GENOMEALIGNMENT_GRAPHCREATOR_H


#include "DataStructures/Read.h"
#include "DataStructures/Graph.h"

class GraphCreator {
public:
    GraphCreator(vector<Read*> *reads, Graph *G);
    virtual ~GraphCreator();

    virtual void startAlignmentGraphCreation() = 0;

    /**
     * Sets alignTo bit at position id to given value
     * @param id
     * @param val
     */
    virtual void setAlignTo( int id, bool val ){ alignTo[id] = val; }
    void setAlignTo( VB & aT ){ alignTo = aT; }

    /**
     * Sets alignFrom bit at position id to given value.
     * @param id
     * @param val
     */
    virtual void setAlignFrom( int id, bool val ){ alignFrom[id] = val; }
    void setAlignFrom( VB & aF ){ alignFrom = aF; }


    bool getAlignFrom(int id){ return alignFrom[id]; }
    bool getAlignTo(int id){ return alignTo[id]; }


    void clear(){}

protected:

//    void addCompRevConnection( Read *r1, Read *r2, int offset ); // if r1 align with r2, then r2_revcomp aligns with r1_revcomp with rightOffset. This function adds that connection
    vector<Read*> *reads;
    Graph * G;
    /**
     * Only reads marked in this vector will be considered as a beginning of an edge in the alignment graph
     */
    VB alignFrom;


    /**
     * Only reads marked in this vector will be considered as an end of an edge in the alignment graph
     */
    VB alignTo;

};


#endif //GENOMEALIGNMENT_GRAPHCREATOR_H
