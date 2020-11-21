//
// Created by sylwester on 12/31/18.
//

#include <GraphCreators/GraphCreatorPrefSuf.h>
#include <Global.h>
#include <Utils/MyUtils.h>
#include <Utils/TimeMeasurer.h>
#include <thread>
#include <AlignmentControllers/AlignmentControllerHybrid.h>
#include <functional>
#include <Utils/WorkloadManager.h>


GraphCreatorPrefSuf::GraphCreatorPrefSuf(vector<Read *> *reads, Graph *G, bool remove_isolated_reads) : GraphCreator(
        reads, G), maxReadLength(0), removeIsolatedReadsBeforeReversingGraph(remove_isolated_reads),
                                                                                                        bitsetChecksCount(
                                                                                                                0),
                                                                                                        bitsetCheckEdgesRemoved(
                                                                                                                0),
                                                                                                        goodPrefSufChecks(
                                                                                                                0),
                                                                                                        prefSufChecks(
                                                                                                                0),
                                                                                                        goodBitsetChecksCount(
                                                                                                                0),
                                                                                                        edgeRemoveAndAddOperations(
                                                                                                                0),
                                                                                                        toRemove(
                                                                                                                Params::THREADS,
                                                                                                                VB(G->size(),
                                                                                                                   false)) {

    calculateMaxReadLength();

    prefixKmers.reserve(G->size());
    prefixKmersAdditional.reserve(G->size());
    suffixKmers.reserve(G->size());
    suffixKmersAdditional.reserve(G->size());

    int lg = (int) log2((double) G->size() / 2);
    prefixKmersBuckets = (1ll << lg);
    if (prefixKmersBuckets * 3ll > G->size()) prefixKmersBuckets >>= 1;
//    prefixKmersBuckets = MyUtils::getNearestLowerPrime((double) G->size() / 3);
    prefixKmersInBuckets = vector<vector<unsigned> >(prefixKmersBuckets);

}

void GraphCreatorPrefSuf::calculateMaxReadLength() {
    for (Read *r : *reads) {
        if (r != nullptr) maxReadLength = max(maxReadLength, r->size());
    }
}

GraphCreatorPrefSuf::~GraphCreatorPrefSuf() {
    clear();
}

void GraphCreatorPrefSuf::clear() {
    vector<unsigned long long>().swap(prefixKmers);
    vector<ADDITIONAL_HASH_TYPE>().swap(prefixKmersAdditional);

    vector<unsigned long long>().swap(suffixKmers);
    vector<ADDITIONAL_HASH_TYPE>().swap(suffixKmersAdditional);

    vector<vector<unsigned> >().swap(prefixKmersInBuckets);

}

void GraphCreatorPrefSuf::startAlignmentGraphCreation() {
    TimeMeasurer::startMeasurement(TimeMeasurer::GRAPH_CREATOR);

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


    cerr << "maxReadLength = " << maxReadLength << endl;
    maxReadLength = min(maxReadLength, 500);

    while (currentPrefSufLength <= maxReadLength) {
        nextPrefSufIteration();
        LL edges = G->countEdges();
        cerr << "After Iteration " << currentPrefSufLength << " / " << maxReadLength <<
             ".  There are already " << edges << " edges in the graph   ->   avg degree "
             << ((double) edges / (double) G->size()) << endl;
    }
    cerr << endl;


    clear();
    if (removeIsolatedReadsBeforeReversingGraph) Global::removeIsolatedReads();

    G->reverseGraphInPlace();

    if (GATHER_STATISTICS) {
        DEBUG(prefSufChecks);
        DEBUG(goodPrefSufChecks);
        DEBUG(bitsetChecksCount);
        DEBUG(goodBitsetChecksCount);
        DEBUG(bitsetCheckEdgesRemoved);
        DEBUG(edgeRemoveAndAddOperations);
        G->writeBasicStatistics();

    }

    Params::MIN_OVERLAP_AREA = oldMOA;
    Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT = oldMOCFA;
    Params::MIN_OVERLAP_RATE = oldMOR;
    Params::MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = oldMOR_ACLER;

    TimeMeasurer::stopMeasurement(TimeMeasurer::GRAPH_CREATOR);
}


void GraphCreatorPrefSuf::createInitialState() {

    unsigned long long dummyKmer = (unsigned long long) (-1); // -1 is just the maximal unsigned long long value
    ADDITIONAL_HASH_TYPE dummyKmerAdditional = (ADDITIONAL_HASH_TYPE) (-1); // -1 is just the maximal unsigned value
    for (unsigned i = 0; i < G->size(); i++) {
        if ((*reads)[i] != nullptr) {
            prefixKmers.emplace_back(0);
            prefixKmersAdditional.emplace_back(0);
            suffixKmers.emplace_back(0);
            suffixKmersAdditional.emplace_back(0);
        } else {
            prefixKmers.push_back(dummyKmer);
            prefixKmersAdditional.push_back(dummyKmerAdditional);
            suffixKmers.push_back(dummyKmer);
            suffixKmersAdditional.push_back(dummyKmerAdditional);
        }
    }

    vector<std::thread> parallelJobs;
    parallelJobs.reserve(Params::THREADS);

    int W = (int) ceil((double) G->size() / Params::THREADS);
    for (int i = 1; i < Params::THREADS; i++) {
        int a = min(i * W, G->size() - 1);
        int b = min((i + 1) * W - 1, G->size() - 1);

        parallelJobs.emplace_back([=] { createInitialStateJob(a, b, i); });
    }

    createInitialStateJob(0, W - 1, 0);


    for (auto &p : parallelJobs) p.join();


    currentPrefSufLength = 0;
    prefHashFactor = 1;
    prefHashFactorAdditional = 1;
    for (int l = 0; l < Params::MIN_OVERLAP_PREF_SUF - 1; l++) {
        currentPrefSufLength++;
        prefHashFactor <<= 2;
        prefHashFactorAdditional <<= 2;

        if (prefHashFactor >= Params::MAX_HASH_CONSIDERED) prefHashFactor %= Params::MAX_HASH_CONSIDERED;
        if (prefHashFactorAdditional >= MAX_ADDITIONAL_HASH) prefHashFactorAdditional %= MAX_ADDITIONAL_HASH;
    }

    cerr << "creating initial state ended" << endl;
}

void GraphCreatorPrefSuf::createInitialStateJob(int a, int b, int thread_id) {
    Params::KMER_HASH_TYPE currentPrefSufLength, prefHashFactor, prefHashFactorAdditional;

    int progressCounter = 0;
    for (int i = a; i <= b; i++) {
        if (!(alignFrom[i] || alignTo[i])) continue;
        currentPrefSufLength = 0;
        prefHashFactor = 1;
        prefHashFactorAdditional = 1;

        for (int l = 0; l < Params::MIN_OVERLAP_PREF_SUF - 1; l++) {
            currentPrefSufLength++;

            bool prefUpdated = alignTo[i] ? updatePrefixHash(i, currentPrefSufLength, prefHashFactor,
                                                             prefHashFactorAdditional) : false;
            if (!prefUpdated) alignTo[i] = false;

            bool suffUpdated = alignFrom[i] ? updateSuffixHash(i, currentPrefSufLength) : false;
            if (!suffUpdated) alignFrom[i] = false;

            prefHashFactor <<= 2;
            if (prefHashFactor >= Params::MAX_HASH_CONSIDERED) prefHashFactor %= Params::MAX_HASH_CONSIDERED;

            prefHashFactorAdditional <<= 2;
            if (prefHashFactorAdditional >= MAX_ADDITIONAL_HASH) prefHashFactorAdditional %= MAX_ADDITIONAL_HASH;
        }

        if (thread_id == 0)
            MyUtils::writeProgress(i - a + 1, b - a + 1, progressCounter, "PrefSuf createInitialState progress", 1);
    }

    if (thread_id == 0) cerr << endl;
}

bool GraphCreatorPrefSuf::updatePrefixHash(int id, int currentPrefSufLength, Params::KMER_HASH_TYPE prefHashFactor,
                                           ADDITIONAL_HASH_TYPE prefHashFactorAdditional) {
    if (currentPrefSufLength > (*reads)[id]->size()) return false;

    prefixKmers[id] += (*(*reads)[id])[currentPrefSufLength - 1] * prefHashFactor;
    if (prefixKmers[id] >= Params::MAX_HASH_CONSIDERED) prefixKmers[id] %= Params::MAX_HASH_CONSIDERED;

    prefixKmersAdditional[id] += (*(*reads)[id])[currentPrefSufLength - 1] * prefHashFactorAdditional;
    if (prefixKmersAdditional[id] >= MAX_ADDITIONAL_HASH) prefixKmersAdditional[id] %= MAX_ADDITIONAL_HASH;
    return true;
}

bool GraphCreatorPrefSuf::updateSuffixHash(int id, int currentPrefSufLength) {
    if (currentPrefSufLength > (*reads)[id]->size() - Params::MIN_OFFSET_FOR_ALIGNMENT) return false;

    suffixKmers[id] <<= 2;
    suffixKmers[id] += (*(*reads)[id])[(*reads)[id]->size() - currentPrefSufLength];
    if (suffixKmers[id] >= Params::MAX_HASH_CONSIDERED) suffixKmers[id] %= Params::MAX_HASH_CONSIDERED;

    suffixKmersAdditional[id] <<= 2;
    suffixKmersAdditional[id] += (*(*reads)[id])[(*reads)[id]->size() - currentPrefSufLength];
    if (suffixKmersAdditional[id] >= MAX_ADDITIONAL_HASH) suffixKmersAdditional[id] %= MAX_ADDITIONAL_HASH;
    return true;
}

void GraphCreatorPrefSuf::nextPrefSufIteration() {
    currentPrefSufLength++;


    vector<std::thread> parallelJobs;
    parallelJobs.reserve(Params::THREADS);

    int W = (int) ceil((double) G->size() / Params::THREADS);
    {
        for (int i = 1; i < Params::THREADS; i++) { // UPDATING PREFIXES
            int a = min(i * W, G->size() - 1);
            int b = min((i + 1) * W - 1, G->size() - 1);
            parallelJobs.emplace_back([=] { updatePrexihHashJob(a, b, i); });
        }
        updatePrexihHashJob(0, W - 1, 0);
        for (auto &p : parallelJobs) p.join();
    }


    {
        W = (int) ceil(prefixKmersBuckets / Params::THREADS);
        parallelJobs.clear();
        for (int i = 1; i < Params::THREADS; i++) { // PLACING KMERS INTO BUCKETS
            int a = min(i * W, prefixKmersBuckets - 1);
            int b = min((i + 1) * W - 1, prefixKmersBuckets - 1);
            parallelJobs.emplace_back([=] { removeKmersFromBucketsJob(a, b, i); });
        }
        removeKmersFromBucketsJob(0, W - 1, 0);
        for (auto &p : parallelJobs) p.join();
    }


    {
        parallelJobs.clear();
        W = (int) ceil((double) G->size() / Params::THREADS);
        for (int i = 1; i < Params::THREADS; i++) { // PLACING KMERS INTO BUCKETS
            int a = min(i * W, G->size() - 1);
            int b = min((i + 1) * W - 1, G->size() - 1);
            parallelJobs.emplace_back([=] { putKmersIntoBucketsJob(a, b, i); });
        }
        putKmersIntoBucketsJob(0, W - 1, 0);
        for (auto &p : parallelJobs) p.join();

        // this sorting should make processor branch prediction better
        /* WorkloadManager::parallelBlockExecution(0, prefixKmersBuckets-1, 3*Params::THREADS, Params::THREADS, [=](unsigned a, unsigned b, unsigned id){
             for(unsigned i=a; i<=b; i++) sort( prefixKmersInBuckets[i].begin(), prefixKmersInBuckets[i].end() );
         });*/
    }


    if (currentPrefSufLength == Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP) {

        G->reverseGraphInPlace();
        G->retainOnlySmallestOffset();
        G->writeBasicStatistics();

        MyUtils::process_mem_usage();
//        Params::THREADS = 1; // #TEST
    }


    {

        int blocks = 50 * Params::THREADS;
        WorkloadManager::parallelBlockExecution(0, G->size() - 1, blocks, Params::THREADS,
                                                [=](unsigned a, unsigned b, unsigned id) {
                                                    nextPrefSufIterationJobAddEdges(a, b, id);
                                                });
    }


    prefHashFactor <<= 2;
    if (prefHashFactor >= Params::MAX_HASH_CONSIDERED) prefHashFactor %= Params::MAX_HASH_CONSIDERED;

    prefHashFactorAdditional <<= 2;
    if (prefHashFactorAdditional >= MAX_ADDITIONAL_HASH) prefHashFactorAdditional %= MAX_ADDITIONAL_HASH;

}

void GraphCreatorPrefSuf::removeKmersFromBucketsJob(int a, int b, int thread_id) {
    for (int i = a; i <= b; i++) {
        vector<unsigned>().swap(prefixKmersInBuckets[i]);
    }
}

void GraphCreatorPrefSuf::putKmersIntoBucketsJob(int a, int b, int thread_id) {
    for (unsigned i = a; i <= b; i++) {
        if (!alignTo[i]) continue;
        int ind = prefixKmers[i] & (prefixKmersBuckets - 1);
//        int ind = prefixKmers[i] % prefixKmersBuckets;
        G->lockNode(ind);
        prefixKmersInBuckets[ind].push_back(i);
        G->unlockNode(ind);
    }
}


void GraphCreatorPrefSuf::writeState() {
    cerr << "currentPrefSufLength: " << currentPrefSufLength << endl;
    cerr << "prefHashFactor = " << prefHashFactor << endl;
    cerr << "maxReadLength = " << maxReadLength << endl;

    cerr << "Prefix kmers:" << endl;
    for (auto a : prefixKmers) cerr << a << endl;

    cerr << endl << "suffixKmers:" << endl;
    for (auto a : suffixKmers) cerr << a << endl;
}

void GraphCreatorPrefSuf::updatePrexihHashJob(int a, int b, int thread_id) {
    for (int i = a; i <= b; i++) {
        if (alignTo[i]) {
            bool prefixUpdated = updatePrefixHash(i, currentPrefSufLength, prefHashFactor, prefHashFactorAdditional);
            if (!prefixUpdated) alignTo[i] = false;
        }
    }
}

void GraphCreatorPrefSuf::nextPrefSufIterationJobAddEdges(int a, int b, int thread_id) {

    vector<unsigned> toR;
    toR.reserve(1000);

    const unsigned MAX_BLOCKS = 10;
    vector<Bitset> temp_bitsets;
    for (int i = 0; i < MAX_BLOCKS; i++) temp_bitsets.emplace_back(Bitset::BLOCK_SIZE * (i + 1));

    auto blNum = [](unsigned offset) { return offset >> Bitset::BLOCK_OFFSET; };
    auto indInBl = [](unsigned index) { return (index) & Bitset::MOD_MODIFIER; };


    for (int suffId = a; suffId <= b; suffId++) {
        bool sufUpdated = false;
        if (alignFrom[suffId]) sufUpdated = updateSuffixHash(suffId, currentPrefSufLength);

        if (sufUpdated) {

            Read *rB = (*reads)[suffId];

            const int b = suffixKmers[suffId] & (prefixKmersBuckets - 1);
            const int offset = (*reads)[suffId]->size() - currentPrefSufLength;

            auto suffHash = suffixKmers[suffId];
            auto suffHashAdditional = suffixKmersAdditional[suffId];

            if (GATHER_STATISTICS) prefSufChecks += prefixKmersInBuckets[b].size();

            for (const unsigned prefId : prefixKmersInBuckets[b]) {
                if (prefId != suffId && prefixKmers[prefId] == suffHash &&
                    prefixKmersAdditional[prefId] == suffHashAdditional) {

                    // uncomment the following section to create full errorless graph
//                    {G->lockNode(prefId); G->pushDirectedEdge(prefId, suffId, offset);  G->unlockNode(prefId); continue;}

                    if (GATHER_STATISTICS) goodPrefSufChecks++;

                    if (Read::calculateReadOverlap((*reads)[suffId], (*reads)[prefId], offset) < currentPrefSufLength)
                        continue; // this line here prohibits included alignment

                    if (currentPrefSufLength < Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP) {

                        // adding normal edges, will have to reverse it later
                        if ((*G)[suffId].size() == SOES) (*G)[suffId].erase((*G)[suffId].begin());
                        G->pushDirectedEdge(suffId, prefId, offset);

                    } else {


                        if (offset > 0) { // this can be checked - maybe i should not consider if offset = 0

                            G->lockNode(prefId);
                            auto neighborhood_list = (*G)[prefId];
                            G->unlockNode(prefId);// #TEST
                            if (GATHER_STATISTICS) bitsetChecksCount += neighborhood_list.size();


                            for (PII &p : neighborhood_list) {
                                const int A = p.first;

                                const int offsetDiff = p.second - offset;


                                if (offsetDiff < 0)
                                    continue; // this is here to prevent checking short edges that were added earlier
                                if (A == suffId)
                                    continue; // if A == B then there will be no edge (A,B) and thus i do not want to remove (A,C) from graph

                                if (GATHER_STATISTICS) goodBitsetChecksCount++;

                                Read *rA = (*reads)[A];

                                if (Read::getRightOffset(rA, rB, offsetDiff) < 0) continue;

                                bool removeEdge;


                                unsigned begBlock = blNum(offsetDiff << 1);
                                unsigned endBlock = blNum(p.second << 1);
                                if (endBlock - begBlock + 1 < MAX_BLOCKS) {
                                    Bitset *temp = &temp_bitsets[endBlock - begBlock];
                                    Bitset *aSeq = &rA->getSequence();
                                    for (int j = begBlock; j <= endBlock; j++)
                                        temp->setBlock(j - begBlock, aSeq->getBlock(j));
                                    (*temp) <<= indInBl(offsetDiff << 1);
                                    Bitset *bSeq = &rB->getSequence();
                                    removeEdge = (!bSeq->mismatchBounded(*temp, offset << 1));
                                } else {
//                                    cerr << "siemka" << endl;
                                    auto bs = rA->getSequence();
                                    bs <<= (offsetDiff << 1);
                                    removeEdge = (!rB->getSequence().mismatchBounded(bs,
                                                                                     (((int) rA->size() - offsetDiff)
                                                                                             << 1)));
                                }

                                if (removeEdge) {
                                    toRemove[thread_id][A] = true;
                                    toR.push_back(A);
                                    if (GATHER_STATISTICS) bitsetCheckEdgesRemoved++;
                                }
                            }
                        }

                        toR.push_back(suffId);
                        toRemove[thread_id][suffId] = true; // this should be with unordered_set<int> version of toRemove

                        G->lockNode(prefId); // #TEST

                        for (unsigned j = (*G)[prefId].size() - 1; j != (unsigned) -1; j--) {
                            if (GATHER_STATISTICS) edgeRemoveAndAddOperations += (*G)[prefId].size();
                            if (toRemove[thread_id][((*G)[prefId][j].first)]) {
                                swap((*G)[prefId][j], (*G)[prefId].back());
                                (*G)[prefId].pop_back();
                            }
                        }

                        for (auto x : toR) toRemove[thread_id][x] = false; // clearing nodes that are to be removed
                        toR.clear();

                        G->pushDirectedEdge(prefId, suffId,
                                            offset); // by adding B to [toRemove], we ensure, that B will not be present in G[C].
                        // If B was present in G[C], then the new connection is better, so it would replace B in the list anyway.
                        if (GATHER_STATISTICS) edgeRemoveAndAddOperations += (*G)[prefId].size();

                        G->unlockNode(prefId);
                    }
                }
            }
        } else alignFrom[suffId] = false;
    }
}




