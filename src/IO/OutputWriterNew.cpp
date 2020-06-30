//
// Created by sylwester on 3/7/19.
//

#include <IO/OutputWriterNew.h>
#include <Utils/MyUtils.h>

#include "IO/OutputWriterNew.h"

OutputWriterNew::OutputWriterNew(Graph * G, vector<Contig *> contigs) {
    this->G = G;
    this->contigs = contigs;
}



void OutputWriterNew::writeContigs() {
    string traverse_type = "_traverse" + to_string( Params::TRAVERSE_TYPE );
    if( Params::TRAVERSE_TYPE == Params::TRAVERSE_GREEDY ) traverse_type += "_ccscl" + to_string( Params::CONTIG_CREATOR_SHORT_CYCLE_LENGTH );

    string rsoe = ( to_string( Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP ) + "-" + to_string( Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN ) );

    Params::outStream.open( Params::outStreamFileName /*+ "_nrpcp" + to_string(Params::NEW_READS_PER_CONTIG_PERCENTAGE) + "_mopp" + to_string(Params::MAX_OFFSET_PARALLEL_PATHS)
        + "_modb" + to_string(Params::MAX_OFFSET_DANGLING_BRANCHES) + rsoe + traverse_type +  ".contigs.fasta"*/ );

    cout.rdbuf( Params::outStream.rdbuf() );

    wasInContig = VB( G->size(), false );

    sort( contigs.begin(), contigs.end(), [](auto &c1, auto &c2){ return c1->size() > c2->size(); } ); // after this contigs should be ordered from longest to shortest.

    int progressCounter = 0;
    int ile = 0;

    int contigsWithEndsInForks = 0;
    int totalContigs = 0;

    cerr << "In OutputWriterNew there are " << contigs.size() << " contigs to write before filtering" << endl;
    for( Contig * ctg : contigs ){


//        cerr << endl << ctg->getId() << endl;


        if( filterContig( ctg ) == false ){
//            cerr  << "Contig " << ctg->getId() << " is rejected" << endl;
            continue;
        }

//        cerr << "Marking contained reads" << endl;
        for( auto p : ctg->getContainedReads() ){
            Read * r = p.first;
//            cerr << "\tr->getId(): " << r->getId() << endl;
            wasInContig[r->getId()] = true;
            if( Params::ADD_COMP_REV_READS ) wasInContig[ r->getIdOfCompRevRead() ] = true;
        }

//        cerr << "Marked contained reads" << endl;

        string s = ctg->getSequenceAsString();
        if( Params::RNA ){
            for( int i=0; i<s.size(); i++ ) if( s[i] == 'T' ) s[i] = 'U';
        }

        if( ctg->endsInFork() ) contigsWithEndsInForks++;
        totalContigs++;

        cout << ">contig_id=" << ctg->getId() << "_length=" << s.size() << endl;

//        writeContig(s,150);
        writeContig(s,1e9);
//        cout << "wasInContig: "; WRITE(wasInContig);

        MyUtils::writeProgress(ile++, contigs.size(), progressCounter, "writing contigs",1 );
    }


    cerr << endl << contigsWithEndsInForks << " out of " << totalContigs << " contigs have ends in a fork" << endl;

}


void OutputWriterNew::writeContig(string &s, int charsPerLine) {
    for( int i=0; i<s.size(); i+= charsPerLine ){
        cout << s.substr( i, charsPerLine ) << endl;
//        cout << s.substr( 50, s.size()-100 ) << endl;
    }
}

bool OutputWriterNew::filterContig(Contig *ctg) {

    if( ctg->size() < Params::CONTIG_MIN_OUTPUT_LENGTH ){
//        cerr << "size of contig: " << ctg->size() << " CONTIG_MIN_OUTPUT_LENGTH = " << Params::CONTIG_MIN_OUTPUT_LENGTH << ", returning false" << endl;
        return false;
    }


    int allReadsInContig = ctg->getContainedReads().size();

    int newReadsInContig = 0;
//    int newReadsInContig = 1;
    for( auto p : ctg->getContainedReads() ){
        Read * r = p.first;
        if( !wasInContig[ r->getId() ] ) newReadsInContig++;
//        wasInContig[r->getId()] = true;
//        wasInContig[ r->getIdOfCompRevRead() ] = true;
    }

    double ratio = (double)newReadsInContig / allReadsInContig;

//    cerr << "considering contig " << *ctg << endl << "  allReadsInContig: " << allReadsInContig << endl << "newReadsInContig: "
//    << newReadsInContig << endl << "ratio: " << ratio << endl;

//    cerr << "New reads in contig: " << newReadsInContig << "   all reads in contig: " << allReadsInContig << endl;

    if( 100*ratio < Params::NEW_READS_PER_CONTIG_PERCENTAGE ){
//        cerr << "ratio = " << ratio << "   returning false" << endl;
        return false;
    }

    return true;
}

VI OutputWriterNew::getContigsLengths() {
    VI res;
    wasInContig = VB( G->size(), false );

    sort( contigs.begin(), contigs.end(), [](auto &c1, auto &c2){ return c1->size() > c2->size(); } ); // after this contigs should be ordered from longest to shortest.

    for( Contig * ctg : contigs ){
        if( filterContig( ctg ) == false ) continue;

        for( auto p : ctg->getContainedReads() ){
            Read * r = p.first;
            wasInContig[r->getId()] = true;
            if( Params::ADD_COMP_REV_READS ) wasInContig[ r->getIdOfCompRevRead() ] = true;
        }

        res.push_back( ctg->size() );
    }

    return res;
}

vector<Contig *> OutputWriterNew::filterContigs() {
    cerr << "Filtering contigs, contigs.size() = " << contigs.size() << endl;

    vector<Contig*> res;

    wasInContig = VB( G->size(), false );

    sort( contigs.begin(), contigs.end(), [](auto &c1, auto &c2){ return c1->size() > c2->size(); } ); // after this contigs should be ordered from longest to shortest.

    int progressCounter = 0;
    int ile = 0;
    int id = 0;

    for( Contig * ctg : contigs ){
        ile++;
        if( filterContig( ctg ) == false ){
//            cerr << "Contig filter returned false" << endl;
            continue;
        }

        for( auto p : ctg->getContainedReads() ){
            Read * r = p.first;
            wasInContig[r->getId()] = true;
            if( Params::ADD_COMP_REV_READS ) wasInContig[ r->getIdOfCompRevRead() ] = true;
        }

        ctg->setId(id++);
        res.push_back(ctg);


        MyUtils::writeProgress(ile, contigs.size(), progressCounter, "filtering contigs",1 );
    }

    return res;

}

void OutputWriterNew::writeContigsNoFilter(vector<Contig *> contigs) {
    Params::outStream.open( Params::outStreamFileName /*+ "_nrpcp" + to_string(Params::NEW_READS_PER_CONTIG_PERCENTAGE) + "_mopp" + to_string(Params::MAX_OFFSET_PARALLEL_PATHS)
        + "_modb" + to_string(Params::MAX_OFFSET_DANGLING_BRANCHES) + rsoe + traverse_type*/ +  ".contigs.fasta" );

    cout.rdbuf( Params::outStream.rdbuf() );

    for( auto t : contigs ){
        string s = t->getSequenceAsString();
        cout << ">contig_id=" << t->getId() << "_length=" << s.size() << endl;
        writeContig( s, 1e9 );
    }
}
