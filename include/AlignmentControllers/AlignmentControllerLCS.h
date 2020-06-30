/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AlignmentControllerLCS.h
 * Author: sylwester
 *
 * Created on November 30, 2018, 1:39 PM
 */

#ifndef ALIGNMENTCONTROLLERLCS_H
#define ALIGNMENTCONTROLLERLCS_H

#include "AlignmentController.h"


class AlignmentControllerLCS : public AlignmentController{
public:
    AlignmentControllerLCS( );
    AlignmentControllerLCS(const AlignmentControllerLCS& orig);
    virtual ~AlignmentControllerLCS();

    
    bool canAlign(Read * r1, Read * r2, int offset) override;
    VVI getFoundLCS(){ return foundLCS; }
private:
    
    VVI lcsdata;
    VVI foundLCS; // foundLCS[0] and foundLCS[1] are a sequence of positions of consequent elements that form LCS. So (*r1)[foundLCS[0][i]] == (*r2)[ foundLCS[1][i] ].
    int calculateLCS(Read* r1, Read* r2, int offset);
};

#endif /* ALIGNMENTCONTROLLERLCS_H */

