//
// Created by sylwester on 12/20/18.
//

#ifndef GENOMEALIGNMENT_GRAPHCREATORKMERBASED_H
#define GENOMEALIGNMENT_GRAPHCREATORKMERBASED_H

#include "GraphCreator.h"
#include "DataStructures/Kmer.h"
#include "Params.h"
#include<cmath>
#include "Utils/MyUtils.h"

#include "StatisticsGenerators/StatisticsGeneratorBigData.h"
#include "AlignmentControllers/AlignmentController.h"
#include "AlignmentControllers/AlignmentControllerHybrid.h"

#include<map>



typedef long long LL;



class GraphCreatorKmerBased : public GraphCreator {

public:

    GraphCreatorKmerBased(vector<Read*> *reads, Graph *G) : GraphCreator(reads,G) { alignmentController = new AlignmentControllerHybrid(); }
    virtual ~GraphCreatorKmerBased();

    /*
     * starts creating alignment graph. There are Params::BUCKETS iterations. In each iterations all kmers with hashes from given bucket are considered, sorted by indexInRead, decreasing.
     * Then for each present hash i take all kmers with that hash and check if they can align.
     */
    virtual void startAlignmentGraphCreation() override;

    static int calculateReadOverlap( Read *r1, Read *r2, int offset ){ return min( r1->size(), r2->size() + offset ) - offset;  }


    virtual void createAlignmentsForKmers( vector<Kmer> & kmers, int p, int q, int thread_id = 0 ) = 0; // given kmers with the same hash, i create graph alignment connections for those kmers that are indexed between p and q inclusive
    void testSerialization();

    virtual GraphCreatorKmerBased * clone() = 0;
protected:

    AlignmentController * alignmentController;

    /*
     * Returns kmers for given bucket (here a bucket is a part of all kmers that will be created. It is used to save space by processing kmers with hashes only from given bucket.
     * The end result will be the same.
     */
    vector<vector<Kmer>> getKmersForBucket(int bucket);

    void getKmersForBucketJob(int a, int b, int thread_id, int bucket, vector<vector<vector<Kmer> > > &bucketKmersInThreads);

    /**
     * Helper function to parallelly invoke createAlignmentForKmers function.
     * In this function there are parallelly created Params::THREADS new GraphCreator's with GraphCreator gc = *this;
     * Each of them parallelly runs createAlignmentForKmers function.
     */
    void createAlignmentForKmersJobNewGC(vector<vector<Kmer>> &kmers, int thread_id);

    /*
     * This can be moved to separate GraphCreatorKmerBased class instance - e.g. GraphCreatorSwats. Is uses fast (non quadratic) algorithm based on offset branches.
     */

    /*
     * Serializes kmers to disk into buckets separate files.
     */
    void serializeKmers( int buckets );

    /*
     * deserializes kmers from given bucket
     */
    vector<Kmer> deserializeKmers(int bucket);

    string getSerializationFilenameForBucket( int bucket );/*{
                                            string s = Params::TEST_NAME + "_bucket" + to_string(bucket+1);
                                            if( !Read::priorities.empty() ){ s += "_prior"; for( auto a : Read::priorities ) s += to_string(a); }
                                            s += ".kmers.bin";  return s; }*/
    void removeSerializationFiles(int buckets);
    bool isAlreadySerialized();

    void sortBucketsJob( vector< vector<Kmer> > & kmers, int a, int b, int thread_id );

private:



    void moveKmersToOneVectorJob(int a, int b, int thread_id, vector< vector< vector<Kmer> > > & kmers);



//    GraphCreatorPairwiseKmer *pairwiseKmerCreator;

    //  map<LL,LL> sameKmersSizes;
 //   StatisticsGeneratorBigData statGen;



};


#endif //GENOMEALIGNMENT_GRAPHCREATORKMERBASED_H
