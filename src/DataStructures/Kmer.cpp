/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Kmer.cpp
 * Author: sylwester
 * 
 * Created on November 23, 2018, 7:59 PM
 */

#include <DataStructures/Kmer.h>

#include "DataStructures/Kmer.h"
#include "DataStructures/Read.h"

/**
 * @param r read from which this kmer comes
 * @param ind index of the beginning of the kmer in given read
 * @param hash this is the hash of the kmer. However i initialize hash in Kmer to hash+KMER_HASH_ADDITION, so that ReadBitMapCreator creates different hashes for
 * multiple occurences of the same SHORT kmer in read. Othwerwise if several short kmers would occur in single read they would get the same hash what is not desirable.
 */
Kmer::Kmer(Read *r, int ind, Params::KMER_HASH_TYPE hash, short length) : read(r), indInRead(
        ind) /*hash(hash + Params::KMER_HASH_ADDITION) */{
    this->hash = hash;
    if (hash < 0) {
        cerr << "hash = " << hash << " < 0 in Kmer contructor!!!" << endl;
        exit(1);
    }
    this->length = length;
}



Kmer::~Kmer() {
    clear();
}

void Kmer::clear() {
}

int Kmer::size() {
    return (int) length;
}


Kmer &Kmer::operator=(const Kmer &oth) {
    read = oth.read;
    indInRead = oth.indInRead;
    hash = oth.hash;
    length = oth.length;
    return *this;
}


bool Kmer::operator<(const Kmer &oth) const {
    if (hash != oth.hash) return hash < oth.hash;
    else if (indInRead != oth.indInRead) return indInRead > oth.indInRead;
    else if (read->size() != oth.read->size()) return read->size() < oth.read->size();
    else return false;

}


string Kmer::getKmerAsString() {
    string s = "";
    for (int i = 0; i < size(); i++) {
        Read *r = read;
        int nukl = (*r)[indInRead + i];
        s += Params::getNuklAsString(nukl);
    }
    return s;
}

ostream &operator<<(ostream &str, Kmer &k) {
    str << "Kmer: " << k.getKmerAsString() << "  with (LL)hash  " << (LL) k.hash << "  at index " << k.indInRead
        << "   with length " << k.length << "  in read " << (*k.read);
    return str;
}

