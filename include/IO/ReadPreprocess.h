//
// Created by sylwester on 11/8/19.
//

#ifndef GENOMEALIGNMENT_READPREPROCESS_H
#define GENOMEALIGNMENT_READPREPROCESS_H

#include "Global.h"


class ReadPreprocess{
public:
    /**
     *
     * @return id's of reads that are a prefix of other reads. If there are two or more exactly the same reads, then the one with greates id will be retained, other will be marked as prefix reads.
     */
    VB getPrefixReads();

    /**
     *
     * @return vector containing pointers to Global::READS, but ordered lexicographically.
     */
    vector<Read*> getSortedReads();

    /**
     *
     * @param reads
     * @return LCP array for given, lexicogarphically sorted reads
     */
    vector<short> getLCP(vector<Read *> &reads);


private:


};

#endif //GENOMEALIGNMENT_READPREPROCESS_H
