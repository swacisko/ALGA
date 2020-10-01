/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AlignmentController.h
 * Author: sylwester
 *
 * Created on November 23, 2018, 8:02 PM
 */

#ifndef ALIGNMENTCONTROLLER_H
#define ALIGNMENTCONTROLLER_H

#include "DataStructures/Read.h"
#include "DataStructures/Kmer.h"
#include "Params.h"
#include<cmath>
#include "Utils/MyUtils.h"
#include "DataStructures/Graph.h"


typedef long long LL;

/*
 Class responsible for creating alignment graph.
 */
class AlignmentController {
public:

    virtual ~AlignmentController() {};

    virtual bool canAlign(Read *r1, Read *r2,
                          int offset) = 0; // returns true, if reads r1 and r2 can be aligned if r2 starts at position offset of r1

//    static int calculateReadOverlap( Read *r1, Read *r2, int offset ){ return min( r1->size(), r2->size() + offset ) - offset;  }
//    static int getRightOffset( Read *r1, Read *r2, int offset ){ return  }
};

#endif /* ALIGNMENTCONTROLLER_H */

