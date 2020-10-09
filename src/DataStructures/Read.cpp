/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Read.cpp
 * Author: sylwester
 * 
 * Created on November 23, 2018, 7:59 PM
 */

#include <DataStructures/Read.h>
#include <numeric>
#include <Utils/MyUtils.h>
#include <Global.h>

#include "DataStructures/Read.h"

Read::Read(int id, string seq) : id(id) {

    createSequence(seq);
}

Read::Read(const Read &orig) : id(orig.id), sequence(orig.sequence) {
}

Read::~Read() {
    clear();
}

int Read::size() {
    return (sequence.size() >> 1);
}

void Read::clear() {
}


void Read::createSequence(string &s) {
    sequence = Bitset((int) (s.size()) << 1);

    for (unsigned i = 0; i < s.size(); i++) {

        switch (s[i]) {
            case 'A': {
                break;
            }
            case 'C': { // C = 01
                sequence.set(i << 1, true);
                break;
            }
            case 'G': { // G = 10
                sequence.set((i << 1) + 1, true);
                break;
            }
            case 'T': { // T = 11
                sequence.set(i << 1, true);
                sequence.set((i << 1) + 1, true);
                break;
            }
            default: { // N = 00
                break;
            }
        }

    }
}

vector<Kmer> Read::getKmers(int length) {
    if (length > size()) return vector<Kmer>();
    if (Params::USE_LI) return getLIKmers(priorities, Params::LI_KMER_LENGTH, Params::LI_KMER_INTERVALS);


    int p = 0;
    int q = 0;
    Params::KMER_HASH_TYPE hash = 0;
    while (q < length) {
        hash <<= 2;
        hash += (*this)[q];

        if (hash >= Params::MAX_HASH_CONSIDERED) hash %= Params::MAX_HASH_CONSIDERED;

        q++;
    }

    Params::KMER_HASH_TYPE factor = 1;
    for (int i = 0; i < length - 1; i++) {
        factor <<= 2; // after this should be 4^(KMER_LENGTH-1)
        if (factor >= Params::MAX_HASH_CONSIDERED) factor %= Params::MAX_HASH_CONSIDERED;
    }


    vector<Kmer> kmers;
    kmers.emplace_back(this, p, hash, length);

    while (q < size()) {
        hash -= factor * (*this)[p];
        if (hash >= Params::MAX_HASH_CONSIDERED || hash < 0) {
            hash %= Params::MAX_HASH_CONSIDERED;
            if (hash < 0) hash += Params::MAX_HASH_CONSIDERED;
        }
        hash <<= 2;
        hash += (*this)[q];
        if (hash >= Params::MAX_HASH_CONSIDERED) hash %= Params::MAX_HASH_CONSIDERED;

        p++;
        q++;
        kmers.emplace_back(this, p, hash, length);
    }


    return kmers;
}


vector<Params::KMER_HASH_TYPE> Read::getKmerHashes(int kmerLength) {
    vector<Kmer> kmers = getKmers(kmerLength);
    vector<Params::KMER_HASH_TYPE> hashes;
    for (Kmer k : kmers) hashes.push_back(k.hash);
    kmers.clear();

    return hashes;
}


Params::NUKL_TYPE Read::operator[](int pos) {
    return ((sequence)[(pos << 1) + 1] << 1) + (sequence)[pos << 1];
}

string Read::getSequenceAsString() {
    string s = "";
    for (int i = 0; i < size(); i++) {
        int nukl = (*this)[i];
        s += Params::getNuklAsString(nukl);
    }
    return s;
}

ostream &operator<<(ostream &str, Read &r) {
    str << r.getSequenceAsString() << "   Read: " << setw(5) << r.id << "   length: " << r.size();
    return str;
}

vector<Kmer> Read::getLIKmers(VI priorities, int length, int intervals) {
    if (length > size()) {
        cerr << "returning empty vector of LIKmers" << endl;
        exit(1);
        return vector<Kmer>();
    }

    typedef __int128 HASH_TYPE;

    int p = 0;
    int q = 0;
    HASH_TYPE hash = 0;
    while (q < length) {
        hash <<= 2;
        hash += priorities.at(Params::getNuklNumber((*this)[q]));
        q++;
    }

    HASH_TYPE factor = 1;
    for (int i = 0; i < length - 1; i++) {
        factor <<= 2; // after this should be 4^(KMER_LENGTH-1)
    }


    vector<HASH_TYPE> minHashInInterval(intervals, factor
            << 2); // i fill with INFINITY (or at least larger than any hash)

    vector<Kmer> minLexInInterval(
            intervals); // minLexInInterval[i] is the LI kmer in i-th interval. i-th interval start at position i*intervalLength and ends in (i+1)*intervalLength



    minHashInInterval[0] = hash;
    minLexInInterval[0] = Kmer(this, p, hash % (HASH_TYPE) Params::MAX_HASH_CONSIDERED, length);

    int intervalLength = (int) ceil(((double) size() - length + 1) / intervals);


    if (intervalLength <= 0) {
        cerr << "intervalLength: " << intervalLength << endl;
        exit(1);
    }

    int interv;

    while (q < size()) {
        hash -= factor * priorities.at(Params::getNuklNumber((*this)[p]));

        hash <<= 2;
        hash += priorities.at(Params::getNuklNumber((*this)[q]));

        p++;
        q++;

        interv = (int) (p / intervalLength);

        if (interv < 0 || interv >= minHashInInterval.size() || interv >= minLexInInterval.size()) {
            cerr << "WRONG BOUND FOR INTERV IN getLIKmers! ! ! ! ! ! ! " << endl;
            exit(1);
        }

        if (hash < minHashInInterval.at(interv)) {
            minHashInInterval.at(interv) = hash;
            minLexInInterval.at(interv) = Kmer(this, p, hash % (HASH_TYPE) Params::MAX_HASH_CONSIDERED, length);
        }
    }

    while (minLexInInterval.size() > interv + 1) minLexInInterval.pop_back();




    for (int j = minLexInInterval.size() - 1; j >= 0; j--) {
        if (minLexInInterval[j].size() == 0) {
            swap(minLexInInterval[j], minLexInInterval.back());
            minLexInInterval.pop_back();
        }
    }

    return minLexInInterval;

}


int Read::getIdOfCompRevRead() {
    if (Params::ADD_COMP_REV_READS) {
        if (id & 1) return id - 1;
        else return id + 1;
    }

    return -1;
}

int Read::getIdOfCompRevRead(int id) {
    if (Params::ADD_COMP_REV_READS) {
        if (id & 1) return id - 1;
        else return id + 1;
    }

    return -1;
}

int Read::getIdOfPairedRead() {

    if (Params::ADD_PAIRED_READS) {
        if (Params::ADD_COMP_REV_READS) {
            if ((id & 3) >= 2) return id - 2;
            else return id + 2;
        } else {
            if (id & 1) return id - 1;
            else return id + 1;
        }
    }
    return -1;
}

int Read::getIdOfPairedRead(int id) {
    if (!Global::pairedReadOffset.empty()) {
        switch (Global::pairedReadOffset[id]) {
            case 0: {
                return id;
            }
            case 1: {
                return id + 2;
            }
            case 2: {
                return id - 2;
            }
            default: {
                return id;
            }
        }
    } else {
        if (Params::ADD_PAIRED_READS) {
            if (Params::ADD_COMP_REV_READS) {
                if ((id & 3) >= 2) return id - 2;
                else return id + 2;
            } else {
                if (id & 1) return id - 1;
                else return id + 1;
            }
        }
        return -1;

    }

}


VI Read::priorities(4, 0);

void Read::set(int pos, Params::NUKL_TYPE val) {
    sequence.set(pos << 1, (bool) (val & 1));
    sequence.set((pos << 1) + 1, (bool) (val & 2));
}

void Read::modifySequence(string &s) {
    createSequence(s);
}
