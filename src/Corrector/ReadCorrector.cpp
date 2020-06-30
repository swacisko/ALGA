//
// Created by sylwester on 2/6/20.
//

#include <Utils/MyUtils.h>
#include <thread>
#include <Corrector/ReadCorrector.h>
#include <Global.h>

#include "Corrector/ReadCorrector.h"

ReadCorrector::ReadCorrector(vector<Read *> &reads, int sLength, int bLength) {
    this->reads = &reads;
    smallLength = sLength;
    bigLength = bLength;
    reversedReads = false;
    mutexes = vector<mutex>(FREQS_SIZE);
}



void ReadCorrector::correct() {
    correctInDirection(false);
    cerr << "reversing direction" << endl;
    correctInDirection(true);
}

void ReadCorrector::correctInDirection(bool reversed) {
    reversedReads = reversed;
    createFrequenciesMap();

    DEBUG(frequencies.size());

    applyCorrection();
    reversedReads = false;
}

void ReadCorrector::addReadDataToMap(Read *r, vector<ReadCorrector::MAP_TYPE> &map) {
    if(r == nullptr) return;
    if( r->size() < smallLength + bigLength ) return;

//    cerr << "Processing read " << r->getSequenceAsString() << endl;

    SMALL_TYPE sH = 0;
    SMALL_TYPE smallPow = (1 << ( 2 * smallLength - 2 )); // 4 ^ (smallLength-1);

    BIG_TYPE bH = 0;
    BIG_TYPE bigPow = ( 1ll << ( 2 * bigLength - 2 ) );

    for( int i=0; i<smallLength; i++ ){
        sH <<= 2;
        sH += accessReadPosition(r,i);
    }
    for( int i=smallLength; i < smallLength + bigLength; i++ ){
        bH <<= 2;
        bH += accessReadPosition(r,i);

        while( bH >= Params::MAX_HASH_CONSIDERED ) bH -= Params::MAX_HASH_CONSIDERED;
    }


    int p = smallLength;
    int q = smallLength + bigLength;

    int bucket = bH & (FREQS_SIZE-1);
    mutexes[bucket].lock();
    map[ bucket ][bH][sH]++;
    mutexes[bucket].unlock();

    while( q < r->size() ){

//        cerr << "\r(p,q) = (" << p << "," << q << ")" << flush;

        sH -= smallPow * accessReadPosition( r, p-smallLength );
        sH <<= 2;
        sH += accessReadPosition(r,p);

        bH -= bigPow * accessReadPosition(r,q-bigLength);
        if( bH < 0 ){
            bH %= Params::MAX_HASH_CONSIDERED;
            if( bH < 0 ) bH += Params::MAX_HASH_CONSIDERED;
        }
        bH <<= 2;
        bH += accessReadPosition(r,q);
        while( bH >= Params::MAX_HASH_CONSIDERED ) bH -= Params::MAX_HASH_CONSIDERED;

        bucket = bH & (FREQS_SIZE-1);
        mutexes[bucket].lock();
        map[ bucket ][bH][sH]++;
        mutexes[bucket].unlock();

        p++;
        q++;
    }

//    cerr << "\tprocessed" << endl;

}


void ReadCorrector::createFrequenciesMap() {
    vector<MAP_TYPE>(FREQS_SIZE).swap(frequencies);

//    vector< MAP_TYPE > freqs( FREQS_SIZE );
    auto helper = [=](int a, int b, int thread_id){
        int progressCounter = 0;
        for( int i=a; i<=b; i++ ){
            if( (*reads)[i] == nullptr ) continue;
            addReadDataToMap((*reads)[i], frequencies);
            if( thread_id == 0 ) MyUtils::writeProgress(i+1, b-a+1, progressCounter, "creating frequencies",1 );
        }
    };


    vector<thread> parallelJobs;
    int W = (int) ceil( (double) reads->size() / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){
        int a = min( i*W, (int)reads->size()-1 );
        int b = min( (i+1)*W-1, (int)reads->size()-1 );
        parallelJobs.push_back( thread( [=] { helper(a, b,i); } ) );
    }
    helper(0, W - 1, 0);
    for( auto & p : parallelJobs ) p.join();



    cerr << "All reads added to map, proceeding to merging maps" << endl;


   /* swap( frequencies, freqs[0]);
    for( int i=1; i<FREQS_SIZE; i++ ){
        for( auto itx : freqs[i] ){
            for( auto ity : itx.second ){
                frequencies[ itx.first ][ ity.first ] += ity.second;
            }
        }

        MAP_TYPE().swap( freqs[i] );
    }*/

//    debugFrequencies();

    cerr << "Retaining only candidates in frequencies, candidateThreshold = " << candidateThreshold << endl;


    auto retainer = [=](int a, int b){
        for( int i=a; i<=b; i++ ){
            auto &freqs = frequencies[i];
            for( auto &itx : freqs ){
                vector<SMALL_TYPE> toRemove;
                for( auto ity : itx.second ){
                    if( ity.second < candidateThreshold ) toRemove.push_back(ity.first);
                }
                for(auto tr : toRemove) itx.second.erase(tr);
            }

            vector<BIG_TYPE> toRemove;
            for(auto itx : freqs) if( itx.second.empty() ) toRemove.push_back(itx.first);
            for(auto tr : toRemove) freqs.erase(tr);

        }
    };

    parallelJobs.clear();
    W = (int) ceil( (double) FREQS_SIZE / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){
        int a = min( i*W, FREQS_SIZE-1 );
        int b = min( (i+1)*W-1, FREQS_SIZE-1 );
        parallelJobs.push_back( thread( [=] { retainer(a, b); } ) );
    }
    retainer(0, W - 1);
    for( auto & p : parallelJobs ) p.join();

    /*for(auto &freqs : frequencies){

        for( auto &itx : freqs ){
            vector<SMALL_TYPE> toRemove;
            for( auto ity : itx.second ){
                if( ity.second < candidateThreshold ) toRemove.push_back(ity.first);
            }
            for(auto tr : toRemove) itx.second.erase(tr);
        }

        vector<BIG_TYPE> toRemove;
        for(auto itx : freqs) if( itx.second.empty() ) toRemove.push_back(itx.first);
        for(auto tr : toRemove) freqs.erase(tr);
    }*/


    cerr << "maps merged" << endl;
//    debugFrequencies();
}



void ReadCorrector::applyCorrection() {

    auto applyCorrectionToReads = [=](int a, int b, int thread_id){
        int progressCounter = 0;
        for(int i=a; i<=b; i++){
            if( (*reads)[i] == nullptr ) continue;
            applyCorrectionToRead( (*reads)[i] );
            if( thread_id == 0 ) MyUtils::writeProgress(i-a+1, b-a+1, progressCounter, "applying correction",1 );
        }
    };

    vector<thread> parallelJobs;
    int W = (int) ceil( (double) reads->size() / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){
        int a = min( i*W, (int)reads->size()-1 );
        int b = min( (i+1)*W-1, (int)reads->size()-1 );

        parallelJobs.push_back( thread( [=] { applyCorrectionToReads(a, b,i); } ) );

        cerr << "(a,b) = (" << a << "," << b << ")" << endl;
    }

    applyCorrectionToReads(0, W - 1, 0);

    for( auto & p : parallelJobs ) p.join();
}

void ReadCorrector::applyCorrectionToRead(Read *r) {
    if(r == nullptr) return;
    if( r->size() < smallLength + bigLength ) return;

//    string rs =  r->getSequenceAsString();
//    if( reversedReads ) reverse(rs.begin(), rs.end());
//    cerr << "Correcting read " << rs << endl;

    SMALL_TYPE sH = 0;
    SMALL_TYPE smallPow = (1 << ( 2 * smallLength - 2 )); // 4 ^ (smallLength-1);

    BIG_TYPE bH = 0;
    BIG_TYPE bigPow = ( 1ll << ( 2 * bigLength - 2 ) );

    for( int i=0; i<smallLength; i++ ){
        sH <<= 2;
        sH += accessReadPosition(r,i);
    }
    for( int i=smallLength; i < smallLength + bigLength; i++ ){
        bH <<= 2;
        bH += accessReadPosition(r,i);

        while( bH >= Params::MAX_HASH_CONSIDERED ) bH -= Params::MAX_HASH_CONSIDERED;
    }


    int p = smallLength;
    int q = smallLength + bigLength;

    auto correctLocal = [=, &r, &p,&q, &sH, &bH](){

//        cerr << "correctLocal: sMer: " << hashToString(sH,smallLength) << "    tMer: " << hashToString(bH, bigLength) << endl;

//        if( frequencies.count(bH) == 0 ) return; // no bigLength kmer - returning (it could have been removed)
        if( frequencies[ bH & (FREQS_SIZE-1) ].count(bH) == 0 ) return; // no bigLength kmer - returning (it could have been removed)

        if( frequencies[ bH & (FREQS_SIZE-1) ][bH].count(sH) ) return; // if sH can be paired with bH in at least candidateThreshold kmers


//        cerr << "Found read to correct" << endl;
//        for(int i=0; i<r->size(); i++) cerr << Params::getNuklAsString( accessReadPosition(r,i) );
//        cerr << endl;
//        for( int i=0; i<p-smallLength; i++ ) cerr << ".";
//        for(int i=p-smallLength; i < q; i++) cerr << Params::getNuklAsString( accessReadPosition(r,i) );
//        for(int i=q; i<r->size(); i++) cerr << ".";
//        cerr << endl;


        SMALL_TYPE closest = -1;
        int minDst = 1e9;

        for( auto ity : frequencies[ bH & (FREQS_SIZE-1) ][bH] ){

            int sMer = ity.first;
            int dst = 0;

            int mask = 3;
            bool sameBoundaries = true;

            for( int i=0; i<smallLength; i++){
                int sMerPos = ( sMer & (mask << (2*i)) ) >> (2*i);
                int readPos = accessReadPosition(r,p-smallLength+i);
                if( sMerPos != readPos ){
                    dst++;

                    // if sMer has different boundaries, then i do not change that - i change only SNPS that are inside the sMer, unless the sMer if at the beginning (or end) of the read)
                    if( (i == 0 || i == smallLength-1) && ( p > smallLength ) ){
//                        cerr << "sH boundaries differ, not correcting" << endl;
//                                cerr << "correctLocal: sMer: " << hashToString(sH,smallLength) << "    tMer: " << hashToString(bH, bigLength) << endl;

                        sameBoundaries = false;
                        break;
                    }
                }
            }

            if( dst < minDst && sameBoundaries ){
                minDst = dst;
                closest = sMer;
            }
        }

        int MAX_SNPS_TO_CORRECT = 1;

        if( minDst > MAX_SNPS_TO_CORRECT ){
//            cerr << "minDst = " << minDst << ",  do not correct" << endl;
            return; /* i do not want to correct more than 2 positions of each sMer*/
        }

//        cerr << "Correcting first " << smallLength << " positions of " << (smallLength+bigLength) << "-mer to: " << hashToString(closest,smallLength) << endl;

        for(int i=0; i<smallLength; i++){
            int mask = 3;
            int sMerPos = ( closest & (mask << (2*i)) ) >> (2*i);
            setReadAtPosition(r,p-1-i,sMerPos);
//            cerr << endl;
        }

//        cerr << "before change sH = " << sH << endl;
        sH = closest;
//        cerr << "After change sH = " << sH << endl;

//        cerr << "Now r:" << endl << *r << endl;


//        exit(1);
    };


    correctLocal();

    while( q < r->size() ){

        sH -= smallPow * accessReadPosition( r, p-smallLength );
        sH <<= 2;
        sH += accessReadPosition(r,p);

        bH -= bigPow * accessReadPosition(r,q-bigLength);
        if( bH < 0 ){
            bH %= Params::MAX_HASH_CONSIDERED;
            if( bH < 0 ) bH += Params::MAX_HASH_CONSIDERED;
        }
        bH <<= 2;
        bH += accessReadPosition(r,q);
        while( bH >= Params::MAX_HASH_CONSIDERED ) bH -= Params::MAX_HASH_CONSIDERED;


        p++;
        q++;

        correctLocal();
    }

}

int ReadCorrector::accessReadPosition(Read *r, int pos) {
    if( reversedReads ){
        return (*r)[r->size()-pos-1];
    }else{
        return (*r)[pos];
    }
}


void ReadCorrector::setReadAtPosition(Read *r, int pos, int val) {
    if( reversedReads ){
        r->set(r->size()-1-pos,val);
    }else{
//        cerr << "Before setting to " << val << ": " << *r << endl;
        r->set(pos,val);
//        cerr << "After setting  to " << val << ": " << *r << endl;
    }
}

void ReadCorrector::debugFrequencies() {

    for(auto freqs : frequencies){
        if( freqs.empty() ) continue;
        for( auto itx : freqs ){
            cerr << hashToString(itx.first, bigLength) << endl;
            for( auto ity : itx.second ){
                cerr << "\t" << hashToString(ity.first, smallLength) << ": " << ity.second << endl;
            }
        }

    }

}

string ReadCorrector::hashToString(LL val, int l) {
    LL mask = 3;
    string res = "";
    for( int i=0; i<l; i++ ){
        res += Params::getNuklAsString( (int)( (val & ( mask << (2*i) ) ) >> 2*i ) );
    }
    reverse(res.begin(), res.end());
    return res;
}



void ReadCorrector::test() {
    Params::THREADS = 3;

    vector<string> str = { "TCTAAAAAAAAAAAAAACCCCCCCCCCCCCCCCTTTTTTTTTTT",
                               "AAAAAAAAAAAAACCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGG",
                                 "AAAAAAAAAAACCCCCCCGCACCCCCCTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGACG",
                                    "GAAAAAAACCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGACGAACTATC",
                                      "AAAAAACCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGACGAACTATCAC",
                                 "AAAAAAAAAAACCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGACGAACTATCACG",
                                        "AAAACCCCCCCCCTTCCCCCTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGACGAACTATCACGG" };

    vector<Read*> reads;
    for(auto s : str){
        reads.push_back( new Read(0,s) );
    }


    cerr << "Before correction" << endl;
    for(auto r : reads){
        cerr << *r << endl;
    }

    cerr << endl;

    ReadCorrector cr(reads,5,15);
    cr.correct();

    cerr << "corrected" << endl;
    for(auto r : reads){
        cerr << *r << endl;
    }

    exit(1);

}


