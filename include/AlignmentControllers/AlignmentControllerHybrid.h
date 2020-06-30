/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AlignmentControllerHybrid.h
 * Author: sylwester
 *
 * Created on November 30, 2018, 1:39 PM
 */

#ifndef ALIGNMENTCONTROLLERHYBRID_H
#define ALIGNMENTCONTROLLERHYBRID_H

#include "AlignmentController.h"
#include "AlignmentControllerLCS.h"
#include "AlignmentControllerLowErrorRate.h"


class AlignmentControllerHybrid : public AlignmentController {
public:
    AlignmentControllerHybrid( );
    AlignmentControllerHybrid(const AlignmentControllerHybrid& orig);
    virtual ~AlignmentControllerHybrid();
    
    bool canAlign(Read* r1, Read* r2, int offset ) override;
    long long totalReadAlignments;
    long long lcsReadAlignments;
    long long bitmapAlignments;
    long long lowErrorAlignmentsApproved;
    long long lowErrorChecks;
    
private:
    AlignmentControllerLCS *lcsController;
    AlignmentControllerLowErrorRate *lerController;



    
};

#endif /* ALIGNMENTCONTROLLERHYBRID_H */

