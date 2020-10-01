//
// Created by sylwester on 12/31/18.
//

#include <GraphCreators/GraphCreatorPrefSuf.h>
#include <Global.h>
#include <Utils/MyUtils.h>
#include <Utils/TimeMeasurer.h>
#include <thread>
#include <AlignmentControllers/AlignmentControllerHybrid.h>
//#include <unordered_map>
#include <functional>


GraphCreatorPrefSuf::GraphCreatorPrefSuf(vector<Read *> *reads, Graph *G) : GraphCreator(reads,G), maxReadLength(0) {
    calculateMaxReadLength();
//    smallOverlapEdges = VVPII(G->size());

    smallOverlapEdges = vector< pair<unsigned,unsigned>* >(G->size());
    for(int i=0; i<G->size(); i++){
        smallOverlapEdges[i] = new pair<unsigned,unsigned>[SOES];
        for(int j=0; j<SOES; j++) smallOverlapEdges[i][j] = {-1,-1};
    }

//    smallOverlapEdges = vector< pair<unsigned,unsigned>[SOES] >(G->size());
//    for(int i=0; i<G->size(); i++) for(int j=0; j<SOES; j++) smallOverlapEdges[i][j] = {-1,-1};

//    for( auto & v : smallOverlapEdges ) v.reserve( Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN+1 );
}

void GraphCreatorPrefSuf::calculateMaxReadLength() {
    for( Read* r : *reads ){
        if( r != nullptr ) maxReadLength = max( maxReadLength, r->size() );
    }
}

GraphCreatorPrefSuf::~GraphCreatorPrefSuf() {
    clear();
}

void GraphCreatorPrefSuf::clear() {
//    for( int i=0; i<prefixKmers.size(); i++ ) if( prefixKmers[i] != nullptr ) delete prefixKmers[i];
    vector<Kmer>().swap( prefixKmers );

//    for( int i=0; i<suffixKmers.size(); i++ ) if( suffixKmers[i] != nullptr ) delete suffixKmers[i];
    vector<Kmer>().swap( suffixKmers );

    vector<vector<Kmer*> >().swap( prefixKmersInBuckets );

}

void GraphCreatorPrefSuf::startAlignmentGraphCreation() {
    TimeMeasurer::startMeasurement( TimeMeasurer::GRAPH_CREATOR );

    int oldMOA = Params::MIN_OVERLAP_AREA = Params::MIN_OVERLAP_PREF_SUF;
    int oldMOCFA = Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT; // here i get INF-1 to avoid accidental specific behaviour for some functions.
    int oldMOR = Params::MIN_OVERLAP_RATE;
    int oldMOR_ACLER = Params::MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR;

    Params::MIN_OVERLAP_AREA = Params::MIN_OVERLAP_PREF_SUF; // exactl value
    Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT = 90; // % of length. Should be enough to guarantee that no positive check will return false in alignment controller hybrid.
    Params::MIN_OVERLAP_RATE = 100; // % of length
    Params::MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = 100; // % of length

    currentPrefSufLength = 0;

    createInitialState();

//    return;
//    exit(1);

    cerr << "maxReadLength = " << maxReadLength << endl;
    maxReadLength = min( maxReadLength, 500 );

    while( currentPrefSufLength <= maxReadLength ){
        nextPrefSufIteration();
        LL edges = G->countEdges();
        cerr << "After Iteration " << currentPrefSufLength << " / " << maxReadLength <<
            ".  There are already " << edges << " edges in the graph   ->   avg degree " << ( (double) edges / (double)G->size() ) << endl;
    }
    cerr << endl;



    clear();
//    if( reads->size() > 1e6 ) Global::removeIsolatedReads();

    G->reverseGraph();




    Params::MIN_OVERLAP_AREA = oldMOA;
    Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT = oldMOCFA;
    Params::MIN_OVERLAP_RATE = oldMOR;
    Params::MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = oldMOR_ACLER;

    TimeMeasurer::stopMeasurement( TimeMeasurer::GRAPH_CREATOR );
}




void GraphCreatorPrefSuf::createInitialState() {
    prefixKmers.reserve( G->size() );

    Kmer dummyKmer(nullptr,0,0,0);
    for( int i=0; i<G->size(); i++ ){
//        if( (*reads)[i] != nullptr ) prefixKmers.push_back( new Kmer( (*reads)[i],0,(Params::KMER_HASH_TYPE )0,0 ) );
//        else prefixKmers.push_back( new Kmer( nullptr, 0, (Params::KMER_HASH_TYPE )0, 0) );

        if( (*reads)[i] != nullptr ) prefixKmers.emplace_back(  (*reads)[i],0,(Params::KMER_HASH_TYPE )0,0 );
        else prefixKmers.push_back( dummyKmer );
    }

    suffixKmers.reserve( G->size() );
    for( int i=0; i<G->size(); i++ ){
//        if( (*reads)[i] != nullptr ) suffixKmers.push_back( new Kmer( (*reads)[i], (*reads)[i]->size(), (Params::KMER_HASH_TYPE )0, 0) );
//        else suffixKmers.push_back( new Kmer( nullptr, 0, (Params::KMER_HASH_TYPE )0, 0) );

        if( (*reads)[i] != nullptr ) suffixKmers.emplace_back( (*reads)[i], (*reads)[i]->size(), (Params::KMER_HASH_TYPE )0, 0 );
        else suffixKmers.push_back( dummyKmer );
    }

//    prefixKmersBuckets = MyUtils::getNearestLowerPrime(G->size() );
    prefixKmersBuckets = MyUtils::getNearestLowerPrime(max(100, G->size() / 3 ) );

    DEBUG(G->size());
    DEBUG(prefixKmersBuckets);


    prefixKmersInBuckets = vector< vector<Kmer*> >( prefixKmersBuckets );


    vector<std::thread> parallelJobs;
    parallelJobs.reserve( Params::THREADS );

    int W = (int) ceil( (double) G->size() / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){
        int a = min( i*W, G->size()-1 );
        int b = min( (i+1)*W-1, G->size()-1 );

//        parallelJobs.push_back( thread( [=] { createInitialStateJob(a, b, i); } ) );
        parallelJobs.emplace_back( [=] { createInitialStateJob(a, b, i); } );
    }

    createInitialStateJob(0, W - 1, 0);


    for( auto & p : parallelJobs ) p.join();



    currentPrefSufLength = 0;
    prefHashFactor = 1;
    for( int l = 0; l < Params::MIN_OVERLAP_PREF_SUF-1; l++ ) {
        currentPrefSufLength++;
        prefHashFactor <<= 2;
        if( prefHashFactor >= Params::MAX_HASH_CONSIDERED ) prefHashFactor %= Params::MAX_HASH_CONSIDERED;
    }

    cerr << "creating initial state ended" << endl;
}

void GraphCreatorPrefSuf::createInitialStateJob(int a, int b, int thread_id) {
    Params::KMER_HASH_TYPE currentPrefSufLength, prefHashFactor;

    int progressCounter = 0;
    for( int i=a; i<=b; i++ ){
        if( !( alignFrom[i] || alignTo[i] ) ) continue;
        currentPrefSufLength = 0;
        prefHashFactor = 1;

        for( int l = 0; l < Params::MIN_OVERLAP_PREF_SUF-1; l++ ){
            currentPrefSufLength++;

            bool prefUpdated = alignTo[i] ? updatePrefixHash(i, currentPrefSufLength, prefHashFactor) : false;
            if( !prefUpdated ) alignTo[i] = false;

            bool suffUpdated = alignFrom[i] ? updateSuffixHash(i, currentPrefSufLength) : false;
            if( !suffUpdated ) alignFrom[i] = false;

            prefHashFactor <<= 2;
            if( prefHashFactor >= Params::MAX_HASH_CONSIDERED ) prefHashFactor %= Params::MAX_HASH_CONSIDERED;

        }

        if( thread_id == 0 ) MyUtils::writeProgress( i-a+1, b-a+1, progressCounter, "PrefSuf createInitialState progress",1 );
    }

    if( thread_id == 0 ) cerr << endl;
}

bool GraphCreatorPrefSuf::updatePrefixHash(int id, int currentPrefSufLength, Params::KMER_HASH_TYPE prefHashFactor) {
    if( currentPrefSufLength > (*reads)[id]->size() ) return false;

    prefixKmers[id].length++;
    prefixKmers[id].hash += (*(*reads)[id])[ currentPrefSufLength-1 ] * prefHashFactor;
    if( prefixKmers[id].hash >= Params::MAX_HASH_CONSIDERED ) prefixKmers[id].hash %= Params::MAX_HASH_CONSIDERED;
    return true;
}

bool GraphCreatorPrefSuf::updateSuffixHash(int id, int currentPrefSufLength) {
    if( currentPrefSufLength > (*reads)[id]->size() - Params::MIN_OFFSET_FOR_ALIGNMENT ) return false;

    suffixKmers[id].length++;
    suffixKmers[id].indInRead--;
    suffixKmers[id].hash <<= 2;
    suffixKmers[id].hash += (*(*reads)[id])[ (*reads)[id]->size() - currentPrefSufLength ];
    if( suffixKmers[id].hash >= Params::MAX_HASH_CONSIDERED ) suffixKmers[id].hash %= Params::MAX_HASH_CONSIDERED;

    return true;
}

void GraphCreatorPrefSuf::nextPrefSufIteration() {
    currentPrefSufLength++;

//    cerr << "Next pref-suf iteration:" << endl;
//    cerr << "Updating prefix hashes" << endl;

//    DEBUG("HERE1");



    vector<std::thread> parallelJobs;
    parallelJobs.reserve(Params::THREADS);

    int W = (int) ceil( (double) G->size() / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){ // UPDATING PREFIXES
        int a = min( i*W, G->size()-1 );
        int b = min( (i+1)*W-1, G->size()-1 );
//        parallelJobs.push_back( thread( [=] { updatePrexihHashJob(a,b,i); } ) );
        parallelJobs.emplace_back( [=] { updatePrexihHashJob(a,b,i); } );
    }
    updatePrexihHashJob(0,W-1,0);
    for( auto & p : parallelJobs ) p.join();


//    cerr << "Removing kmers from buckets" << endl;
    W = (int) ceil( prefixKmersBuckets / Params::THREADS );
    parallelJobs.clear();
    for( int i=1; i<Params::THREADS; i++ ){ // PLACING KMERS INTO BUCKETS
        int a = min( i*W, prefixKmersBuckets-1 );
        int b = min( (i+1)*W-1, prefixKmersBuckets-1 );
//        parallelJobs.push_back( thread( [=] { removeKmersFromBucketsJob(a,b,i); } ) );
        parallelJobs.emplace_back(  [=]{ removeKmersFromBucketsJob(a,b,i);} );
    }
    removeKmersFromBucketsJob(0,W-1,0);
    for( auto & p : parallelJobs ) p.join();


//    cerr << "Putting kmers into buckets" << endl;
    parallelJobs.clear();
    W = (int) ceil( (double) G->size() / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){ // PLACING KMERS INTO BUCKETS
        int a = min( i*W, G->size()-1 );
        int b = min( (i+1)*W-1, G->size()-1 );
//        parallelJobs.push_back( thread( [=] { putKmersIntoBucketsJob(a,b,i); } ) );
        parallelJobs.emplace_back(  [=] { putKmersIntoBucketsJob(a,b,i); } );
    }
    putKmersIntoBucketsJob(0,W-1,0);
    for( auto & p : parallelJobs ) p.join();






    if( currentPrefSufLength == Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP ){
        cerr << "moving small overlap edges to graph" << endl;
        parallelJobs.clear();
        for( int i=1; i<Params::THREADS; i++ ){ // PLACING KMERS INTO BUCKETS
            int a = min( i*W, G->size()-1 );
            int b = min( (i+1)*W-1, G->size()-1 );
//            parallelJobs.push_back( thread( [=] { moveSmallOverlapEdgesToGraphJob(a,b,i); } ) );
            parallelJobs.emplace_back(  [=] { moveSmallOverlapEdgesToGraphJob(a,b,i); } );
        }
        moveSmallOverlapEdgesToGraphJob(0,W-1,0);
        for( auto & p : parallelJobs ) p.join();

//        cerr << "Reversing graph!! CAUTION, this should be used only with adding non-reversed edges in moveSmallOverlapEdgesToGraphJob" << endl;
//        G->reverseGraph();

        G->retainOnlySmallestOffset();

        cerr << "After moving small overlap edges to graph, G has " << G->countEdges() << " edges" << endl;

        MyUtils::process_mem_usage();

        vector<SOES_TYPE>().swap( smallOverlapEdges );
    }



//    cerr << "Adding edges for iteration" << endl;
    parallelJobs.clear();
    for( int i=1; i<Params::THREADS; i++ ){ // PLACING KMERS INTO BUCKETS
        int a = min( i*W, G->size()-1 );
        int b = min( (i+1)*W-1, G->size()-1 );
//        parallelJobs.push_back( thread( [=] { nextPrefSufIterationJobAddEdges(a,b,i); } ) );
        parallelJobs.emplace_back( [=] { nextPrefSufIterationJobAddEdges(a,b,i); } );
    }
    nextPrefSufIterationJobAddEdges(0,W-1,0);
    for( auto & p : parallelJobs ) p.join();


    prefHashFactor <<= 2;
    if( prefHashFactor >= Params::MAX_HASH_CONSIDERED ) prefHashFactor %= Params::MAX_HASH_CONSIDERED;

}

void GraphCreatorPrefSuf::removeKmersFromBucketsJob(int a, int b, int thread_id) {
    for( int i=a; i<=b; i++ ) {
        vector<Kmer*>().swap( prefixKmersInBuckets[i] );
    }
}

void GraphCreatorPrefSuf::putKmersIntoBucketsJob(int a, int b, int thread_id) {
    for( int i=a; i<=b; i++ ){
        if( !alignTo[i] ) continue;
        int ind =  prefixKmers[i].hash % prefixKmersBuckets;
        G->lockNode(ind);
        prefixKmersInBuckets[ ind ].push_back( &prefixKmers[i] );
        G->unlockNode(ind);
    }
}


void GraphCreatorPrefSuf::writeState() {
    cerr << "currentPrefSufLength: " << currentPrefSufLength << endl;
    cerr << "prefHashFactor = " << prefHashFactor << endl;
    cerr << "maxReadLength = " << maxReadLength << endl;

    cerr << "Prefix kmers:" << endl;
    for( auto a : prefixKmers ) cerr << a << endl;

    cerr << endl << "suffixKmers:" << endl;
    for(auto a : suffixKmers) cerr << a << endl;
}

void GraphCreatorPrefSuf::updatePrexihHashJob(int a, int b, int thread_id) {
    for( int i=a; i<=b; i++ ){
        if( alignTo[i] ){
            bool prefixUpdated = updatePrefixHash(i, currentPrefSufLength, prefHashFactor);
            if( !prefixUpdated ) alignTo[i] = false;
        }
    }
}

void GraphCreatorPrefSuf::nextPrefSufIterationJobAddEdges(int a, int b, int thread_id) {
    VPII toRemove;

//    static long long kmersChecked = 0;
//    static long long kmersHashEquals = 0;

    for( int i=a; i<=b; i++ ){
        bool sufUpdated = false;
        if( alignFrom[i] ) sufUpdated = updateSuffixHash(i, currentPrefSufLength);



        if( sufUpdated ){
            Kmer * suff = &suffixKmers[i];
            const int suffId = suff->read->getId(); // should be suffId == i

            const int b = suff->hash % prefixKmersBuckets;

            const int offset = suff->read->size() - currentPrefSufLength;


            for( Kmer * pref : prefixKmersInBuckets[b] ){

//                kmersChecked++;

                int prefId = pref->read->getId();
                if( pref->hash == suff->hash && prefId != suffId  ){

//                    kmersHashEquals++;

                    if( Read::calculateReadOverlap( suff->read, pref->read, offset ) < currentPrefSufLength ) continue; // this line here prohibits included alignment

                        if( currentPrefSufLength < Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP ){

                            int xx = SOES; for( int j=0; j<SOES; j++ ) if( smallOverlapEdges[suffId][j] == pair<unsigned,unsigned>(-1,-1) ){ xx = j; break; }
                            if( xx == SOES ){
                                for(int j=0; j<SOES-1; j++){
                                    smallOverlapEdges[suffId][j] = smallOverlapEdges[suffId][j+1];
                                }
                                xx = SOES-1;
                            }

                            smallOverlapEdges[ suffId ][xx] = {prefId,offset} ; // i add normal edges here


//                            smallOverlapEdges[ suffId ].push_back( {prefId,offset} ); // i add normal edges here
//                            if( smallOverlapEdges[suffId].size() > Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN ){
//                                smallOverlapEdges[suffId].erase( smallOverlapEdges[suffId].begin() );
//
//                                if( smallOverlapEdges[suffId].size() > Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN ){
//                                    cerr << "ERROR in smallOverlapEdges.size() in PS0" << endl; exit(1);
//                                }
//                            }
                        }else{

                            const int C = prefId;
                            const int B = suffId;
                            G->lockNode(C);

                            auto neighborhood_list = (*G)[C];// #TEST
                            G->unlockNode(C);// #TEST

                            if( offset > 0) { // this can be checked - maybe i should not consider if offset = 0

//                                for (PII &p : (*G)[C]) {
                                for (PII &p : neighborhood_list) { // #TEST
                                    const int A = p.first;

                                    const int offsetDiff = p.second - offset;

                                    if( offsetDiff < 0 ) continue; // this is here to prevent checking short edges that were added earlier
                                    if( A == B ) continue; // if A == B then there will be no edge (A,B) and thus i do not want to remove (A,C) from graph

                                    auto bs = (*reads)[A]->getSequence();
                                    bs <<= (offsetDiff<<1);
                                    bool removeEdge = /* offsetDiff > 0 &&*/  Read::getRightOffset( (*reads)[A], (*reads)[B], offsetDiff) >= 0
//                                                                                && ach->canAlign((*reads)[A], (*reads)[B], offsetDiff); // the constraint offsetDiff > 0 may perhaps be avoided
//                                            &&  ( (*reads)[B]->getSequence().mismatch( (*reads)[A]->getSequence() << (offsetDiff<<1) ) >=  (int)(*reads)[A]->size() - offsetDiff   );
                                            // this line above is the same condition as in ach->canAlign(), but done faster for no-error overlap
                                            &&  ( (*reads)[B]->getSequence().mismatch( bs ) >=  (int)(*reads)[A]->size() - offsetDiff   );



                                    if ( removeEdge )  {
                                        toRemove.push_back( {C, A} );
                                    }
                                }
                            }

                            G->lockNode(C); // #TEST
                            for( auto &p : toRemove ){
                                G->removeDirectedEdge( p.first, p.second );
                            }
                            toRemove.clear();
                            G->addDirectedEdge( C,B, offset );

                            G->unlockNode(C);

                        }

                }
            }
        }else alignFrom[i] = false;

    }

//    DEBUG(kmersChecked);
//    DEBUG(kmersHashEquals);

}


void GraphCreatorPrefSuf::moveSmallOverlapEdgesToGraphJob(int a, int b, int thread_id) {
    for( int i=a; i<=b; i++ ){
        Kmer * suff = &suffixKmers[i];
        if( suff == nullptr || suff->read == nullptr ) continue;

        int suffId = suff->read->getId();
        if( currentPrefSufLength == Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP ){
//            for( auto & p : smallOverlapEdges[suffId] ){
            for( int j=0; j<SOES; j++ ){
                auto p = smallOverlapEdges[suffId][j];
                if(p == pair<unsigned,unsigned>(-1, -1) ) break;
                G->lockNode( p.first );
                (*G)[p.first].push_back( { suffId, p.second } ); // i add reverse edges to the graph!!
                G->unlockNode(p.first);
            }
//            pair<unsigned,unsigned>[3].swap( smallOverlapEdges[suffId] );
            if( smallOverlapEdges[suffId] != nullptr ){
                delete[] smallOverlapEdges[suffId];
                smallOverlapEdges[suffId] = nullptr;
            }
        }
    }
}



