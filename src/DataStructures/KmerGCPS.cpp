//
// Created by sylwester on 10/1/20.
//

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   KmerGCPS.cpp
 * Author: sylwester
 *
 * Created on November 23, 2018, 7:59 PM
 */

#include <DataStructures/KmerGCPS.h>
#include "DataStructures/Read.h"

/**
 * @param r read from which this kmer comes
 * @param ind index of the beginning of the kmer in given read
 * @param hash this is the hash of the kmer. However i initialize hash in Kmer to hash+KMER_HASH_ADDITION, so that ReadBitMapCreator creates different hashes for
 * multiple occurences of the same SHORT kmer in read. Othwerwise if several short kmers would occur in single read they would get the same hash what is not desirable.
 */
KmerGCPS::KmerGCPS(unsigned id, unsigned hash) : read_id(id) {
    this->hash = hash;
    if (hash < 0) {
        cerr << "hash = " << hash << " < 0 in Kmer contructor!!!" << endl;
        exit(1);
    }
}


KmerGCPS &KmerGCPS::operator=(const KmerGCPS &oth) {
    read_id = oth.read_id;
    hash = oth.hash;
    return *this;
}


bool KmerGCPS::operator<(const KmerGCPS &oth) const {
    return hash < oth.hash;
}


ostream &operator<<(ostream &str, KmerGCPS &k) {
    str << "Kmer with read_id " << k.read_id << " and hash  " << (LL) k.hash;
    return str;
}

