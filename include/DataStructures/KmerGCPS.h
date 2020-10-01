//
// Created by sylwester on 10/1/20.
//

#ifndef ALGA_KMERGCPS_H
#define ALGA_KMERGCPS_H

#include<iostream>

using namespace std;

class Read;

/**
 * This is just the kmer class that
 */
class KmerGCPS {
public:
    KmerGCPS() : read_id(0), hash(0) {}

    KmerGCPS(unsigned id,
             unsigned hash); // creates Kmer for given read. This Kmer starts at the position ind and has given length

//    Read *read; // pointer to the read this Kmer comes from
    unsigned read_id; // id of the read this kmer comes from
    unsigned long long hash;

    bool operator<(const KmerGCPS &oth) const;

    KmerGCPS &operator=(const KmerGCPS &oth);


    friend ostream &operator<<(ostream &str, KmerGCPS &k);

private:


};

#endif //ALGA_KMERGCPS_H
