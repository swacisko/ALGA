//
// Created by sylwester on 3/4/19.
//

#ifndef GENOMEALIGNMENT_CONTIG_H
#define GENOMEALIGNMENT_CONTIG_H


#include <DataStructures/Read.h>

class Contig : public Read {
public:
    Contig(int id, string s, vector<pair<Read *, int>> &containedReads);

    ~Contig();

    vector<pair<Read *, int> > &getContainedReads();

    void setContainedReads(const vector<pair<Read *, int>> &containedReads);

//    int getIdOfCompRevRead(){ return -1; };
    int getIdOfPairedRead() { return -1; };

    /**
     * This function is used to correct SNPs that can occur when creating string brutally. Here we consider all reads that create this contig and
     * for each posiion take the nucleotide hat occurs in most reads (we take most frequent nucleotide on given position).
     * When this functon is invoked, the @sequence parameter of this contig will be changed to proper one.
     * @return correct string representation of this contig
     */
    string correctSnipsInContig();

    void writeContainedReads();

    /**
     *
     * @return true if this contig ends in a fork, false otherwise
     */
    bool endsInFork() { return ends_in_fork; }

    void setEndsInFork(bool val) { ends_in_fork = val; }

    static int ID_COUNT;

    static void test();

    /**
     * Modifies this contigs sequence to given sequence
     * @param s
     */
    void modifySequence(string &s);


private:

    /**
     * True if this contig ends in fork, false if this contig has end not in a fork
     */
    bool ends_in_fork = false;

    /**
     * containedReads[i] contains i-th read contained in this contig and offset from (i-1)-th read to this read.
     */
    vector<pair<Read *, int> > containedReads;


};


#endif //GENOMEALIGNMENT_CONTIG_H
