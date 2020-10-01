/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AlignmentControllerHybrid.cpp
 * Author: sylwester
 * 
 * Created on November 30, 2018, 1:39 PM
 */

#include "AlignmentControllers/AlignmentControllerHybrid.h"

AlignmentControllerHybrid::AlignmentControllerHybrid() {
    lcsController = new AlignmentControllerLCS();
    lerController = new AlignmentControllerLowErrorRate();

    totalReadAlignments = 0;
    lcsReadAlignments = 0;
    bitmapAlignments = 0;
    lowErrorAlignmentsApproved = 0;
    lowErrorChecks = 0;
}

AlignmentControllerHybrid::AlignmentControllerHybrid(const AlignmentControllerHybrid &orig) : AlignmentController(
        orig) {
    lcsController = new AlignmentControllerLCS();
    lerController = new AlignmentControllerLowErrorRate();

    totalReadAlignments = 0;
    lcsReadAlignments = 0;
    bitmapAlignments = 0;
    lowErrorAlignmentsApproved = 0;
    lowErrorChecks = 0;
}

AlignmentControllerHybrid::~AlignmentControllerHybrid() {
    delete lcsController;
    lcsController = 0;
    delete lerController;
    lerController = 0;
}

bool AlignmentControllerHybrid::canAlign(Read *r1, Read *r2, int offset) {

    try {
        if (100 * offset >
            Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT * r1->size() /*min( r1->size(), r2->size() )*/  )
            return false;

        if (offset < Params::MIN_OFFSET_FOR_ALIGNMENT) return false;

        int overlap = Read::calculateReadOverlap(r1, r2, offset);
        if (overlap < Params::MIN_OVERLAP_AREA) return false;

        if (Read::getRightOffset(r1, r2, offset) < 0) return false;


        totalReadAlignments++;


        if (Params::USE_LCS_LOW_ERROR_FILTER) {
            lowErrorChecks++;


            if (lerController->canAlign(r1, r2, offset)) {
                lowErrorAlignmentsApproved++;
                return true;
            } else if (Params::USE_ACLER_INSTEAD_OF_ACLCS) return false;
        }

        lcsReadAlignments++;
        return lcsController->canAlign(r1, r2, offset);

    }
    catch (int e) {
        cerr << "exception caught in ACH" << endl;
        exit(1);
    }

}



