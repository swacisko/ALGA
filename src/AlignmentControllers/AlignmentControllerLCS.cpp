/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AlignmentControllerLCS.cpp
 * Author: sylwester
 * 
 * Created on November 30, 2018, 1:39 PM
 */

#include <Utils/TimeMeasurer.h>
#include <AlignmentControllers/AlignmentControllerLCS.h>
#include <unordered_map>

#include "AlignmentControllers/AlignmentControllerLCS.h"

AlignmentControllerLCS::AlignmentControllerLCS() {
    // lcsdata = VVI( 200, VI(200,0) );
}

AlignmentControllerLCS::AlignmentControllerLCS(const AlignmentControllerLCS &orig) : AlignmentController(orig) {
}

AlignmentControllerLCS::~AlignmentControllerLCS() {
}

bool AlignmentControllerLCS::canAlign(Read *r1, Read *r2, int offset) {
//    cerr << "can align LCS" << endl;
    if (100 * offset >
        Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT * r1->size() /*min( r1->size(), r2->size() )*/  )
        return false;
//    TimeMeasurer::startMeasurement(TimeMeasurer::ALIGNMENT_CONTROLLER_CAN_ALIGN_LCS);

    int overlap_area = Read::calculateReadOverlap(r1, r2, offset);
    if (overlap_area < Params::MIN_OVERLAP_AREA) return false;

    int overlap = calculateLCS(r1, r2, offset);

    // BELOW LINES ARE TO ENSURE THAT ENDS OF OVERLAP ARE THE SAME. THIS HOWEVER IS DUBIOUS APPLICATION IN LCS, BECAUSE
    // USING LCS WE WANT TO DETECT INDELS, SO ENDS WILL PROBABLY NOT BE THE SAME, EVEN IF OVERLAPS ARE ALMOST THE SAME (SAME BUT WITH INDEL)
//    for( int i=0; i<Params::ALIGNMENT_CONTROLLER_SAME_ENDS_LENGTH; i++ ){
//        if( (*r1)[offset+i] != (*r2)[i] ) return false;
//        if( (*r1)[offset+overlap-1-i] != (*r2)[overlap-1-i] ) return false;
//    }



//    int readOverlap = Read::calculateReadOverlap(r1,r2,offset);
    int readOverlap = r1->calculateReadOverlap(r1, r2, offset);


//    TimeMeasurer::stopMeasurement(TimeMeasurer::ALIGNMENT_CONTROLLER_CAN_ALIGN_LCS);

    if (100 * overlap > Params::MINIMAL_OVERLAP_RATE_FOR_LCS * readOverlap)return true;
    else return false;
}

int AlignmentControllerLCS::calculateLCS(Read *r1, Read *r2, int offset) {
    int M = max(r1->size(), r2->size());


//    if( M > lcsdata.size() ){
//        VVI().swap(lcsdata);
//        lcsdata = VVI( 2*M+5, VI(2*M+5,0) );
//    }

    vector<unordered_map<int, int> > lcsdata(2 * M + 5);


    int E = Params::MAX_ERROR_RATE_FOR_LCS;

    int pBeg = max(0, offset - E);

    for (int p = pBeg; p < r1->size(); p++) {

        int qBeg = max(0, p - offset - E);
        int qEnd = min(r2->size() - 1, p - offset + E);

        if (p == pBeg && p > 0) {
            for (int q = max(0, qBeg - 1); q <= qEnd; q++) {
                lcsdata[p - 1][q] = 0;
                //   if( q > 0 ) lcsdata[p-1][q-1] = 0;
            }
        }
        for (int q = max(0, qBeg - 1); q <= min(r2->size() - 1, qEnd + 1); q++) {
            lcsdata[p][q] = 0; // here i clear the data - it may contain some calculations from previous calculateLCS function invokes.
            //  if( q > 0 ) lcsdata[p][q-1] = 0;
        }

        for (int q = qBeg; q <= qEnd; q++) {
            if ((*r1)[p] == (*r2)[q]) {
                if (p > 0 && q > 0) lcsdata[p][q] = lcsdata[p - 1][q - 1] + 1;
                else lcsdata[p][q] = 1;
            } else {
                if (p > 0 && lcsdata[p][q] < lcsdata[p - 1][q]) lcsdata[p][q] = lcsdata[p - 1][q];
                if (q > 0 && lcsdata[p][q] < lcsdata[p][q - 1]) lcsdata[p][q] = lcsdata[p][q - 1];
            }

        }
    }


//    int p = r1->size()-1;
    int p = min(r1->size() - 1, r2->size() - 1 + offset);

//    int qBeg = maxVal( 0, p-offset-E );
    //   int qEnd = minVal( r2->size()-1, p-offset+E );
    //   int q = minVal( r2->size()-1, qEnd+1 );
    int q = min(r2->size() - 1, p - offset + E);




    //   return lcsdata[p][q]; // this can be done to return only the length of the lcs.



    //*** SECTION BELOW IS TO EXTRACT EXACT POSITION OF LCS
//    foundLCS.clear();
//    foundLCS = VVI(2);
//
//    // this loop will recreate the lcs based on the data in lcsdata.
//    while( p >= pBeg && q >= 0 ){
//        if( (*r1)[p] == (*r2)[q] ){
//            foundLCS[0].push_back(p);
//            foundLCS[1].push_back(q);
//        //    cout << "\tadding, p = " << p << "   q = " << q << "  r1[p] = " << (*r1)[p] << "   r2[q] = " << (*r2)[q] << endl;
//            p--;
//            q--;
//        }else{
//            if( p > pBeg && lcsdata[p][q] == lcsdata[p-1][q] ) p--;
//            else if( q > 0 && lcsdata[p][q] == lcsdata[p][q-1] ) q--;
//            else p--; // here can be q--, no matter.
//        }
//    }
//
//    reverse( foundLCS[0].begin() , foundLCS[0].end() );
//    reverse( foundLCS[1].begin() , foundLCS[1].end() );




    p = min(r1->size() - 1, r2->size() - 1 + offset);
    q = min(r2->size() - 1, p - offset + E);

    return lcsdata[p][q];
}
