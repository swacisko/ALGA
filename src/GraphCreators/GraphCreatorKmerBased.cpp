//
// Created by sylwester on 12/20/18.
//

#include "GraphCreators/GraphCreatorKmerBased.h"


#include <numeric>
#include "GraphSimplifiers/GraphSimplifier.h"
#include "AlignmentControllers/AlignmentController.h"

#include "Global.h"
#include "Utils/TimeMeasurer.h"
#include "StatisticsGenerators/GenomeStatisticsCollector.h"
#include<algorithm>
#include <GraphCreators/GraphCreatorKmerBased.h>
#include <fstream>
#include <thread>
//


GraphCreatorKmerBased::~GraphCreatorKmerBased() {
    delete alignmentController;
    alignmentController = 0;
}


void GraphCreatorKmerBased::startAlignmentGraphCreation() {
    TimeMeasurer::startMeasurement( TimeMeasurer::GRAPH_CREATOR );

    int BUCKETS = 1;
    for( int i=0; i<BUCKETS; i++ ) {

        string s = "Bucket #" + MyUtils::toString<int>(i + 1);
        cerr << s << flush;
        TimeMeasurer::startMeasurement(s);


        vector<vector<Kmer> > kmers = getKmersForBucket(i);

        long long elements = 0;
        for (auto &k : kmers) elements += k.size();
        cerr << endl << "There are " << elements << " kmers in the bucket." << endl;


        cerr << "Largest bucket has " << max_element( kmers.begin(), kmers.end(),[](auto & v, auto & u){return v.size() < u.size();} )->size() << " elements" << endl;

        cerr << endl << "Proceeding to sorting buckets" << endl;

        TimeMeasurer::startMeasurement(TimeMeasurer::KMER_BUCKETS_SORTING);

        vector<std::thread> parallelJobs;


        int W = (int) ceil( (double) kmers.size() / Params::THREADS);
        for( int i=1; i<Params::THREADS; i++ ){
            int a = i*W;
            int b = min( (i+1)*W-1, (int)kmers.size()-1 );
            parallelJobs.push_back( thread( [=, &kmers] { sortBucketsJob(kmers,a,b,i); } ) );
        }
        sortBucketsJob(kmers,0,W-1,0);
        for( auto & p : parallelJobs ) p.join();
        cerr << endl;

        parallelJobs.clear();

        TimeMeasurer::stopMeasurement(TimeMeasurer::KMER_BUCKETS_SORTING);

        cerr << "Buckets sorted" << endl;

        for (int i = 1; i < Params::THREADS; i++) {
            parallelJobs.push_back(thread([=, &kmers] { createAlignmentForKmersJobNewGC(kmers, i); }));
        }

        createAlignmentForKmersJobNewGC(kmers, 0);
        for (auto &p : parallelJobs) p.join();
        TimeMeasurer::stopMeasurement(s);

        cerr << endl << "Alignments for kmers done" << endl;


    }

    G->retainOnlySmallestOffset(); // this should not make any difference in theory, but i still add it anyway (in case it does in practice)

    TimeMeasurer::stopMeasurement( TimeMeasurer::GRAPH_CREATOR );



}

void GraphCreatorKmerBased::sortBucketsJob(vector< vector<Kmer> > & kmers, int a, int b, int thread_id){


    int progressCounter=0;
    for (int i = a; i <= b; i++) {
        if( kmers[i].size() > 0 ) sort( kmers[i].begin(), kmers[i].end() );

        if( thread_id==0 ) MyUtils::writeProgress(i + 1-a, b-a+1, progressCounter, "sorting buckets progress", 1);

    }

}

void GraphCreatorKmerBased::createAlignmentForKmersJobNewGC(vector<vector<Kmer>> &kmers, int thread_id) {
    GraphCreatorKmerBased *gc = (GraphCreatorKmerBased*) this->clone();

    int p,q;
    int progressCounter = 0;
    int totalSize = 0;
    for( int i=0; i<kmers.size(); i++ ) if( i % Params::THREADS == thread_id ) totalSize += kmers[i].size();
    int ile = 0;

    for( int i=0; i<kmers.size(); i++ ){
        if( i % Params::THREADS != thread_id ) continue;
        p = q = 0;

        while( p < (kmers)[i].size() ) {
            while (q < (kmers)[i].size() && (kmers)[i][q].hash == (kmers)[i][p].hash) q++;

            gc->createAlignmentsForKmers( kmers[i],p,q-1, thread_id );

            ile += q-p;
            p = q;
        }

        if( thread_id == 0 ) MyUtils::writeProgress( ile,totalSize, progressCounter, "bucket processing", 10 );

    }

    delete gc;
    gc = 0;
}



vector<vector<Kmer>> GraphCreatorKmerBased::getKmersForBucket(int bucket) {
    long long BUCKETS_SORT = 1048576ll;  // this is the number of buckets for fast sort of all elements.

    vector< vector< vector<Kmer> > > bucketKmersInThreads( Params::THREADS, vector<vector<Kmer>>(BUCKETS_SORT) );

    vector<std::thread> parallelJobs;

    int W = (int) ceil( (double)reads->size() / Params::THREADS);
    for( int i=1; i<Params::THREADS; i++ ){
        int a = i*W;
        int b = min( (i+1)*W-1, (int)reads->size()-1 );
        parallelJobs.push_back( thread( [=, &bucketKmersInThreads] { getKmersForBucketJob(a,b,i,bucket, bucketKmersInThreads); } ) );
    }

    getKmersForBucketJob(0,W-1,0,bucket,bucketKmersInThreads);
    for( auto & p : parallelJobs ) p.join();


    cerr << endl << "Proceeding to parallel adding kmers to one vector" << endl;


    //*********** PARALLEL ADDING KMERS TO ONE VECTOR
    parallelJobs.clear();
    W = (int) ceil( (double)BUCKETS_SORT / Params::THREADS);
    for( int i=1; i<Params::THREADS; i++ ){
        int a = i*W;
        int b = min( (i+1)*W-1, (int)BUCKETS_SORT-1 );

        parallelJobs.push_back( thread( [=, &bucketKmersInThreads] { moveKmersToOneVectorJob(a,b,i, bucketKmersInThreads); } ) );
    }

    moveKmersToOneVectorJob(0,W-1,0, bucketKmersInThreads);
    for( auto & p : parallelJobs ) p.join();


    return bucketKmersInThreads[0];

}


void GraphCreatorKmerBased::moveKmersToOneVectorJob(int a, int b, int thread_id, vector<vector< vector<Kmer> > > &bucketKmersInThreads) {
    int progressCounter = 0;

    for( int i=a; i<=b; i++ ){
        int totalSize = 0;
        for( int k=1; k<Params::THREADS; k++ ) totalSize += bucketKmersInThreads[k][i].size();
        bucketKmersInThreads[0][i].reserve( bucketKmersInThreads[0][i].size() + totalSize ); // that will be the number of elements after moving

        for( int k=1; k<Params::THREADS; k++ ) {
            bucketKmersInThreads[0][i].insert( bucketKmersInThreads[0][i].end(), bucketKmersInThreads[k][i].begin(), bucketKmersInThreads[k][i].end() );
            vector<Kmer>().swap( bucketKmersInThreads[k][i] );
        }

        if( thread_id == 0 ) MyUtils::writeProgress(i+1,b-a+1, progressCounter, "moving kmers to one vector",1);
    }
}

void GraphCreatorKmerBased::getKmersForBucketJob(int a, int b, int thread_id, int bucket, vector<vector<vector<Kmer> > > &bucketKmersInThreads) {
    Params::KMER_HASH_TYPE MAX_HASH = Params::MAX_HASH_CONSIDERED;

    Params::KMER_HASH_TYPE blockSize = (MAX_HASH) + 1;

    Params::KMER_HASH_TYPE A = blockSize * bucket;
    Params::KMER_HASH_TYPE B = blockSize*(bucket+1)-1;

    long long BUCKETS_SORT = 1048576ll;  // 1048576ll; // this is the number of buckets for fast sort of all elements.

    int readsProcessed = 0;
    int progressCounter = 0;
    if( thread_id == 0 ) cerr << endl << "creating kmers for bucket" << flush;

    for( int i=a; i<=b; i++ ){
        Read *r = (*reads)[i];
        if( r == nullptr ) continue;
        vector<Kmer> readKmers;

        if( !( alignFrom.at(r->getId()) || alignTo.at(r->getId()) ) ) continue;


        if( readKmers.empty() ){
            readKmers = r->getKmers(Params::KMER_LENGTH_BUCKET);
        }

        readsProcessed++;

        for( Kmer &k : (readKmers) ){
            if( k.hash >= A && k.hash <= B ){
                int ind = (int)(  (BUCKETS_SORT-1) * ( (double)(k.hash-A)  / (double)(B-A) )  );
                if( ind < 0 || ind >= bucketKmersInThreads[thread_id].size() ){
                    cerr << "ERROR in getKmersForBucket, bucket = " << bucket << "   (LL)A = " << (long long)A << "   (LL)B = " << (long long)B << "  bucketsForSort = " <<
                         BUCKETS_SORT << "  (LL)k.hash = " << (long long)k.hash << "  ind = " << ind << endl;
                    exit(1);
                }
                bucketKmersInThreads[thread_id][ind].push_back(k);

                if(k.size() == 0 ){
                    cerr << "k.size() == 0 in getKmersForBucketJob" << endl;
                    cerr << k << endl;
                    exit(1);
                }
            }
        }

        if( thread_id == 0 ) MyUtils::writeProgress( readsProcessed,Global::READS.size() / Params::THREADS, progressCounter, "creating kmers for bucket",1 );


        vector<Kmer>().swap(readKmers);
    }

    if(thread_id == 0) cerr << endl;
}




