/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Kmer.h
 * Author: sylwester
 *
 * Created on November 23, 2018, 7:59 PM
 */

#ifndef KMER_H
#define KMER_H

#include "Params.h"
#include<iostream>
using namespace std;

class Read;

class Kmer {
public:
    Kmer() : read(nullptr), indInRead(0), hash(0), length(0) {}
    Kmer(Read *r, int ind, Params::KMER_HASH_TYPE hash, short length); // creates Kmer for given read. This Kmer starts at the position ind and has given length
//    Kmer(const Kmer& orig);
    /*virtual*/ ~Kmer();
    
   // void operator=( const Kmer &oth );
    
    void clear();
    int size();


    Read* read; // pointer to the read this Kmer comes from
    Params::KMER_HASH_TYPE hash;
    int indInRead; //Indeks of the beginning of this Kmer in the read.
    short length;


    bool operator<( const Kmer & oth ) const;
    Kmer& operator=( const Kmer & oth );

//    bool operator==( const Kmer & oth ) const;

    string getKmerAsString();

    friend ostream& operator<<( ostream& str, Kmer & k );
private:


};



#endif /* KMER_H */
