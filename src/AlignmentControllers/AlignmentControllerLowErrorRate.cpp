//
// Created by sylwester on 12/16/18.
//

#include <AlignmentControllers/AlignmentControllerLowErrorRate.h>
#include <Global.h>
#include <Utils/TimeMeasurer.h>


AlignmentControllerLowErrorRate::AlignmentControllerLowErrorRate() {
    b1 = b2 = Bitset(1);
}


bool AlignmentControllerLowErrorRate::canAlign(Read *r1, Read *r2, int offset) {
    if (100 * offset >
        Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT * r1->size() /*min( r1->size(), r2->size() ) */ )
        return false;

    int overlap_area = Read::calculateReadOverlap(r1, r2, offset);
    if (overlap_area < Params::MIN_OVERLAP_AREA) return false;

    if (offset < 0) {
        cerr << "offset < 0 in ACLER!" << endl;
        exit(1);
    }

    b1 = r1->getSequence();
    b2 = r2->getSequence();

//    int readOverlap = Read::calculateReadOverlap(r1,r2,offset);
    int readOverlap = r1->calculateReadOverlap(r1, r2, offset);
    b1 <<= (2 * offset);
    b1 ^= b2;
    int sequenceOverlap = (readOverlap << 1) - b1.count(0, (readOverlap << 1) - 1);
    sequenceOverlap >>= 1;


    /**
     * First and last @SAME_ENDS_LENGTH nucleotides of the overlap must be the same.
     */
    int SAME_ENDS_LENGTH = Params::ALIGNMENT_CONTROLLER_SAME_ENDS_LENGTH;
    if (b1.count(0, SAME_ENDS_LENGTH << 1) != 0 ||
        b1.count((readOverlap - SAME_ENDS_LENGTH) << 1, (readOverlap << 1) - 1) != 0)
        return false;

    if (100 * sequenceOverlap >= Params::MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR * readOverlap) return true;
    else return false;
}

void AlignmentControllerLowErrorRate::test() {

    string s1 = "ACTGACTGGGACTGACTTTAAAA";
    string s2 = "GGGGGACTGACTGGGACTGACTTT";
    int offset = 0;

    Read *r1 = new Read(1, s1);
    Read *r2 = new Read(2, s2);

    AlignmentControllerLowErrorRate acler;

    bool can = acler.canAlign(r1, r2, offset);

    cerr << "canAlign: " << can << endl;
}

