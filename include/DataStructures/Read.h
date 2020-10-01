/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Read.h
 * Author: sylwester
 *
 * Created on November 23, 2018, 7:59 PM
 */

#ifndef READ_H
#define READ_H

#include "Bitset.h"
#include "Params.h"
#include "Kmer.h"
#include<string>
#include<iomanip>


class Kmer;

class Read {
public:
    Read(int id, string sequence);

    Read(const Read &orig);

//protected:
//    virtual ~Read();
    ~Read();

public:

    void clear();

    int size();

    Params::NUKL_TYPE operator[](int pos);

    void set(int pos, Params::NUKL_TYPE val);

    void changeId(
            int newId) { id = newId; } // changes id to newID. This is used to remove duplicate reads from all READS and assign new id's to those reads

    bool operator<(Read &oth) { /*if( *sequence != *oth.sequence ) return (*sequence) < (*oth.sequence); else */return
                id < oth.id;
    }

    bool operator==(Read &oth) { /*return (*sequence) == (*oth.sequence);*/ return id == oth.id; }

    Bitset &getSequence() { return sequence; }


    vector<Params::KMER_HASH_TYPE> getKmerHashes(int kmerLength);

    vector<Kmer> getKmers(int length);

    int getId() { return id; }

    void setId(int val) { id = val; }

    string getSequenceAsString();

    friend ostream &operator<<(ostream &str, Read &r);
//    bool operator==( const Read & oth );

    static VI priorities;

    static void writeReadsWithOffset(Read *r1, Read *r2, int offset) {
        cerr << *r1 << endl;
        for (int i = 0; i < offset; i++) cerr << " ";
        cerr << *r2 << endl;
    }

    int
    getIdOfCompRevRead(); // returns id of a read that is reverse complementary of this one or -1 if we did not add reverse complemetary.
    int
    getIdOfPairedRead(); // return id of a read that is paired with this one or -1 if we do not have such information.

    static int getIdOfPairedRead(int id);

    static int getIdOfCompRevRead(
            int id); // returns id of a read that is reverse complementary of this one or -1 if we did not add reverse complemetary.

    /*
     * returns overlap for two given reads.
     */
    static int calculateReadOverlap(Read *r1, Read *r2, int offset) {
        return min(r1->size(), r2->size() + offset) - offset;
    }

    /*
     * returns right offset for given two reads and offset. It is used in adding revcomp connections.
     */
    static int getRightOffset(Read *r1, Read *r2, int offset) { return r2->size() + offset - r1->size(); }


    /**
     * Modifies this contigs sequence to given sequence
     * @param s
     */
    void modifySequence(string &s);


protected:
    Bitset sequence;
    int id;

    void createSequence(string &s);


    vector<Kmer> getLIKmers(VI priorities, int length, int intervals);


};

#endif /* READ_H */

