//
// Created by sylwester on 12/31/18.
//

#include <GraphCreators/GraphCreatorPrefSuf.h>
#include <Global.h>
#include <Utils/MyUtils.h>
#include <Utils/TimeMeasurer.h>
#include <thread>
#include <AlignmentControllers/AlignmentControllerHybrid.h>
#include <unordered_map>
#include <functional>


GraphCreatorPrefSuf::GraphCreatorPrefSuf(vector<Read *> *reads, Graph *G) : GraphCreator(reads,G), maxReadLength(0) {
    calculateMaxReadLength();
//    smallOverlapEdges = VVPII(G->size());
    smallOverlapEdges = vector< pair<unsigned,unsigned>[SOES] >(G->size());
    for(int i=0; i<G->size(); i++) for(int j=0; j<SOES; j++) smallOverlapEdges[i][j] = {-1,-1};
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
        cerr << "Iteration " << currentPrefSufLength << " / " << maxReadLength << ".  There are already " << edges << " edges in the graph   ->   avg degree " << ( (double) edges / (double)G->size() ) << endl;
    }
    cerr << endl;



    clear();
    if( reads->size() > 1e6 ) Global::removeIsolatedReads();

    G->reverseGraph();

//    removeSmallOverlapEdges();



    Params::MIN_OVERLAP_AREA = oldMOA;
    Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT = oldMOCFA;
    Params::MIN_OVERLAP_RATE = oldMOR;
    Params::MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = oldMOR_ACLER;

    TimeMeasurer::stopMeasurement( TimeMeasurer::GRAPH_CREATOR );
}

void GraphCreatorPrefSuf::removeSmallOverlapEdges() {
    cerr << endl << "Removing small overlap edges, G has " << G->countEdges() << " edges" << endl;

    int threshold = Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN;

    function< void(int,int,int) > threadJob = [=](int a, int b, int thread_id){
        for( int i=a; i<=b; i++ ){

            if( (*reads)[i] == nullptr || (*G)[i].empty() ) continue;

            Read* r1 = (*reads)[i];

            sort( (*G)[i].begin(), (*G)[i].end(), [=,&r1]( PII a, PII b ){
                int overlapA = Read::calculateReadOverlap( r1, (*reads)[a.first], a.second );
                int overlapB = Read::calculateReadOverlap( r1, (*reads)[b.first], b.second );

                return overlapA > overlapB;
            } );

            PII back = (*G)[i].back();

            if( (*reads)[ back.first ] == nullptr ){
                cerr << "r2 == nullptr" << endl;
                exit(1);
            }

            Read* r2 = (*reads)[ back.first ];
            int offset = back.second;
            int overlap = Read::calculateReadOverlap( r1,r2,offset );

            if( outdegOverThreshold[i] >= threshold ){

                while( overlap < Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP && !(*G)[i].empty() ){
                    (*G)[i].pop_back();
                    back = (*G)[i].back();
                    r2 = (*reads)[ back.first ];
                    offset = back.second;
                    overlap = Read::calculateReadOverlap( r1,r2,offset );
                }
            }else{
                int neigh = outdegOverThreshold[i];
                int p=0;
                while( p < (*G)[i].size()
                    && Read::calculateReadOverlap( r1, (*reads)[ (*G)[i][p].first ], (*G)[i][p].second ) >= Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP  ) p++;

                while( neigh < Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN && (*G)[i].size() > p ){
                    (*G)[i].pop_back();
                    neigh++;
                }
            }
        }
    };

    vector<thread> parallelJobs;
    int W = (int) ceil( (double) G->size() / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){
        int a = min( i*W, G->size()-1 );
        int b = min( (i+1)*W-1, G->size()-1 );

        parallelJobs.push_back( thread( [=] { threadJob(a, b, i); } ) );
    }

    threadJob(0, W - 1, 0);

    for( auto & p : parallelJobs ) p.join();

    cerr << "Removed, now G has " << G->countEdges() << " edges" << endl;
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

//    prefixKmersBuckets = G->size();
    prefixKmersBuckets = MyUtils::getNearestLowerPrime(G->size());

    DEBUG(G->size());
    DEBUG(prefixKmersBuckets);


    prefixKmersInBuckets = vector< vector<Kmer*> >( prefixKmersBuckets );


    vector<std::thread> parallelJobs;

    int W = (int) ceil( (double) G->size() / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){
        int a = min( i*W, G->size()-1 );
        int b = min( (i+1)*W-1, G->size()-1 );

        parallelJobs.push_back( thread( [=] { createInitialStateJob(a, b, i); } ) );
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

    int W = (int) ceil( (double) G->size() / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){ // UPDATING PREFIXES
        int a = min( i*W, G->size()-1 );
        int b = min( (i+1)*W-1, G->size()-1 );
        parallelJobs.push_back( thread( [=] { updatePrexihHashJob(a,b,i); } ) );
    }
    updatePrexihHashJob(0,W-1,0);
    for( auto & p : parallelJobs ) p.join();


//    cerr << "Removing kmers from buckets" << endl;
    W = (int) ceil( prefixKmersBuckets / Params::THREADS );
    parallelJobs.clear();
    for( int i=1; i<Params::THREADS; i++ ){ // PLACING KMERS INTO BUCKETS
        int a = min( i*W, prefixKmersBuckets-1 );
        int b = min( (i+1)*W-1, prefixKmersBuckets-1 );
        parallelJobs.push_back( thread( [=] { removeKmersFromBucketsJob(a,b,i); } ) );
    }
    removeKmersFromBucketsJob(0,W-1,0);
    for( auto & p : parallelJobs ) p.join();


//    cerr << "Putting kmers into buckets" << endl;
    parallelJobs.clear();
    W = (int) ceil( (double) G->size() / Params::THREADS );
    for( int i=1; i<Params::THREADS; i++ ){ // PLACING KMERS INTO BUCKETS
        int a = min( i*W, G->size()-1 );
        int b = min( (i+1)*W-1, G->size()-1 );
        parallelJobs.push_back( thread( [=] { putKmersIntoBucketsJob(a,b,i); } ) );
    }
    putKmersIntoBucketsJob(0,W-1,0);
    for( auto & p : parallelJobs ) p.join();






    if( currentPrefSufLength == Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP ){
        cerr << "moving small overlap edges to graph" << endl;
        parallelJobs.clear();
        for( int i=1; i<Params::THREADS; i++ ){ // PLACING KMERS INTO BUCKETS
            int a = min( i*W, G->size()-1 );
            int b = min( (i+1)*W-1, G->size()-1 );
            parallelJobs.push_back( thread( [=] { moveSmallOverlapEdgesToGraphJob(a,b,i); } ) );
        }
        moveSmallOverlapEdgesToGraphJob(0,W-1,0);
        for( auto & p : parallelJobs ) p.join();


        G->retainOnlySmallestOffset();

        cerr << "After moving small overlap edges to graph, G has " << G->countEdges() << " edges" << endl;

        MyUtils::process_mem_usage();
        vector<pair<unsigned,unsigned>[SOES]>().swap( smallOverlapEdges );
    }



//    cerr << "Adding edges for iteration" << endl;
    parallelJobs.clear();
    for( int i=1; i<Params::THREADS; i++ ){ // PLACING KMERS INTO BUCKETS
        int a = min( i*W, G->size()-1 );
        int b = min( (i+1)*W-1, G->size()-1 );
        parallelJobs.push_back( thread( [=] { nextPrefSufIterationJobAddEdges(a,b,i); } ) );
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

    for( int i=a; i<=b; i++ ){
        bool sufUpdated = false;
        if( alignFrom[i] ) sufUpdated = updateSuffixHash(i, currentPrefSufLength);



        if( sufUpdated ){
            Kmer * suff = &suffixKmers[i];
            int suffId = suff->read->getId(); // should be suffId == i

            int b = suff->hash % prefixKmersBuckets;

            int offset = suff->read->size() - currentPrefSufLength;


            for( Kmer * pref : prefixKmersInBuckets[b] ){

                int prefId = pref->read->getId();
                if( pref->hash == suff->hash && prefId != suffId  ){

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

                            int C = prefId;
                            int B = suffId;
                            G->lockNode(C);

                            if( offset > 0 ) { // this can be checked - maybe i should not consider if offset = 0

                                for (PII &p : (*G)[C]) {
                                    int A = p.first;

                                    int offsetDiff = p.second - offset;

                                    if( offsetDiff < 0 ) continue; // this is here to prevent checking short edges that were added earlier
                                    if( A == B ) continue; // if A == B then there will be no edge (A,B) and thus i do not want to remove (A,C) from graph

                                    auto bs = (*reads)[A]->getSequence();
                                    bs <<= (offsetDiff<<1);
                                    bool removeEdge = /* offsetDiff > 0 &&*/  Read::getRightOffset( (*reads)[A], (*reads)[B], offsetDiff) >= 0
//                                                                                && ach->canAlign((*reads)[A], (*reads)[B], offsetDiff); // the constraint offsetDiff > 0 may perhaps be avoided
//                                            &&  ( (*reads)[B]->getSequence().mismatch( (*reads)[A]->getSequence() << (offsetDiff<<1) ) >=  (int)(*reads)[A]->size() - offsetDiff   );
                                            &&  ( (*reads)[B]->getSequence().mismatch( bs ) >=  (int)(*reads)[A]->size() - offsetDiff   );
                                            // this line above is the same condition as in ach->canAlign(), but done faster for no-error overlap
;
//                                            bool removeEdge = false;
                                    if ( removeEdge )  {
                                        toRemove.push_back( {C, A} );
                                    }
                                }
                            }

                            for( auto &p : toRemove ){
                                G->removeDirectedEdge( p.first, p.second );
                            }
                            toRemove.clear();
                            G->addDirectedEdge( C,B, offset );

                            G->unlockNode(C);


//                            G->lockNode(B);
//                            if( outdegOverThreshold[B] < 250 ) outdegOverThreshold[B]++;
//                            G->unlockNode(B);
                        }

                }
            }
        }else alignFrom[i] = false;

    }

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
        }
    }
}



