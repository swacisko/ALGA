//
// Created by sylwester on 12/16/18.
//

#ifndef GENOMEALIGNMENT_ALIGNMENTCONTROLLERLOWERRORRATE_H
#define GENOMEALIGNMENT_ALIGNMENTCONTROLLERLOWERRORRATE_H


#include "AlignmentController.h"

/*
 * This alignment controller should work fine for reads with almost none changes (it takes bit & operation on read sequences to check if they align.
 */
class AlignmentControllerLowErrorRate : public AlignmentController {

public:
    AlignmentControllerLowErrorRate();

    bool canAlign(Read *r1, Read *r2, int offset) override;

    static void test();
private:

    Bitset b1;
    Bitset b2;

};


#endif //GENOMEALIGNMENT_ALIGNMENTCONTROLLERLOWERRORRATE_H
