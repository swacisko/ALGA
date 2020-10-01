//
// Created by sylwester on 3/4/19.
//

#include <DataStructures/Contig.h>
#include <Global.h>

#include "DataStructures/Contig.h"

Contig::Contig(int id, string s, vector<pair<Read *, int>> &containedReads) : Read(id, s) {
    this->containedReads = containedReads;
}

Contig::~Contig() {

}


void Contig::setContainedReads(const vector<pair<Read *, int>> &containedReads) {
    Contig::containedReads = containedReads;
}

int Contig::ID_COUNT = 0;

vector<pair<Read *, int> > &Contig::getContainedReads() {
    return containedReads;
}

void Contig::modifySequence(string &s) {
    createSequence(s);
}


string Contig::correctSnipsInContig() {

    string s = "";
    vector<pair<Read *, int> > correctors;
    correctors.emplace_back(containedReads[0].first, 0);
    VI mostFrequent(4);

    containedReads.emplace_back(new Read(-1, s),
                                containedReads.back().first->size()); // here i add any read as a sentry to the end of contained reads. it will be removed later.

    vector<short> freqs;
    for (int i = 1; i < containedReads.size(); i++) {
        int offset = containedReads[i].second;


//        cerr << "next read: " << *containedReads[i].first << endl;

        while (offset > 0) {
            offset--;

            fill(mostFrequent.begin(), mostFrequent.end(), 0);
            for (int k = correctors.size() - 1; k >= 0; k--) {

                Read *r = correctors[k].first;
                int ind = correctors[k].second;

                if (ind >= r->size()) {
                    swap(correctors[k], correctors.back());
                    correctors.pop_back();
                    continue;
                }

                correctors[k].second++;
                mostFrequent[(*r)[ind]]++;
            }

            auto it = max_element(mostFrequent.begin(), mostFrequent.end());
            freqs.push_back(*it);

//            cerr << "\tappending " << Params::getNuklAsString( it - mostFrequent.begin() ) << " to s" << endl;
            s += Params::getNuklAsString(it - mostFrequent.begin());

        }

        if (i < containedReads.size() - 1) correctors.emplace_back(containedReads[i].first, 0);

    }

    containedReads.pop_back();


    int THR = 3;
    int p = 0, q = freqs.size() - 1;
    while (p <= q && freqs[p] <= THR) p++;
    while (p <= q && freqs[q] <= THR) q--;

    s = s.substr(p, q - p + 1);



    /*int TRIM_END_LEFT = 14;
    int TRIM_END_RIGHT = 8;

    if( containedReads.size() > 2 ){
        TRIM_END_LEFT = containedReads[1].second *//*+ containedReads[2].second*//*; // i trim the whole part of the first read

        Read* r1 = containedReads[ getContainedReads().size()-2 ].first;
        Read* r2 = containedReads.back().first;
        int offset = containedReads.back().second;

        TRIM_END_RIGHT = Read::getRightOffset( r1,r2,offset );
    }


    if( s.size() > TRIM_END_LEFT + TRIM_END_RIGHT  ) s = s.substr( TRIM_END_LEFT, s.size() - TRIM_END_LEFT - TRIM_END_RIGHT);
    else s = "A"; // for some reason i do not want empty contig :)*/


//    exit(1);

    createSequence(s);
    return s;
}

void Contig::writeContainedReads() {
    int offset = 0;
    for (int i = 0; i < containedReads.size(); i++) {
        if (i > 0) {
            offset += containedReads[i].second;
            for (int k = 0; k < offset; k++) cerr << " ";
        }
        cerr << *containedReads[i].first << endl;
    }
}

void Contig::test() {

    string A = "ACGGGGTGTGTGTGTAACC";
    string B = "GGGAGTGTGTGTAACCGGTCC";
    string C = "GAGTGTGTGTAACCCGTCC";
    string D = "GAGTGTGTGTAACCGGTCCAAA";
//    string res =  "GGGAGTGTGTGTAACCGGTCC";

    string s = "AAAAAAAAAAAAAAAAAAa";

    vector<pair<Read *, int> > containedReads;
    int id = 0;
    containedReads.push_back({new Read(id++, A), 0});
    containedReads.push_back({new Read(id++, B), 3});
    containedReads.push_back({new Read(id++, C), 2});
    containedReads.push_back({new Read(id++, D), 0});

    Contig *ctg = new Contig(id++, s, containedReads);

    cerr << "before correcting SNPs:" << endl << *ctg << endl;
    ctg->correctSnipsInContig();
    cerr << "after correcting SNPs:" << endl << *ctg << endl;

    exit(1);
}


