//
// Created by sylwester on 3/7/19.
//

#ifndef GENOMEALIGNMENT_OUTPUTWRITERNEW_H
#define GENOMEALIGNMENT_OUTPUTWRITERNEW_H


#include <DataStructures/Contig.h>
#include <DataStructures/Graph.h>

class OutputWriterNew {

public:

    OutputWriterNew(Graph *G, vector<Contig *> contigs);

    void writeContigsNoFilter(vector<Contig *> contigs);

    void writeContigs();

    void writeContig(string &s, int charsPerLine);

    /**
     * Function responsible for filtering contigs.
     * If this function returns true, then given contig will be written to output. Otherwise a contig will not be written to output.
     */
    bool filterContig(Contig *ctg);

    VI getContigsLengths();

    /**
     * Function filters and returns only those contigs that are to be written to output.
     * @return
     */
    vector<Contig *> filterContigs();

private:
    Graph *G;
    vector<Contig *> contigs;

    /**
     * This variable stores information about reads that were already processed in contigs. This is used to decrease duplication ratio.
     */
    VB wasInContig;

};


#endif //GENOMEALIGNMENT_OUTPUTWRITERNEW_H
