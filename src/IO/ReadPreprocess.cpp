//
// Created by sylwester on 11/8/19.
//


#include <IO/ReadPreprocess.h>
#include <future>

#include "IO/ReadPreprocess.h"
#include "Utils/TimeMeasurer.h"

VB ReadPreprocess::getPrefixReads() {

    vector<Read *> reads = getSortedReads();
    vector<VB> markedToRemove(Params::THREADS, VB(Global::READS.size(), false));



    // memory efficient version - uses bitmask instead of returning ids of nodes
    cerr << "Removing prefix reads - calculating LCP" << endl;
    TimeMeasurer::startMeasurement("LCP");
    int N = reads.size();

    auto lcpFun = [&reads, &markedToRemove](int c, int d, int thread_id) {
        for (int i = c; i <= d; i++) {
            int a = i;
            int b = i + 1;
            Read *r1 = reads[a];
            Read *r2 = reads[b];

            int m = min(r1->size(), r2->size());
            int l = 0;
            int BEG = 0;
            if (r1->getSequence().getBlock(0) == r2->getSequence().getBlock(0)) l = BEG = Bitset::BLOCK_SIZE / 2;

            for (int j = BEG; j < m; j++) {
                if ((*r1)[j] == (*r2)[j]) l++;
                else break;
            }

            if (Params::REMOVE_PREF_READS_TYPE == Params::PREF_READS_ONLY_DUPLICATES) {
                if (l == reads[i]->size() && reads[i]->size() == reads[i + 1]->size()) {
//                    markedToRemove[ reads[i]->getId() ] = true;
                    markedToRemove[thread_id][reads[i]->getId()] = true;
                }
            } else if (Params::REMOVE_PREF_READS_TYPE == Params::PREF_READS_ALL_PREFIX_READS && l == reads[i]->size()) {
//                markedToRemove[ reads[i]->getId() ] = true;
                markedToRemove[thread_id][reads[i]->getId()] = true;

                if (reads[i]->size() < reads[i + 1]->size()) {
//                    markedToRemove[ reads[i]->getIdOfCompRevRead() ] = true; // comprev read is a proper suffix of another read, we remove it as well.
                    markedToRemove[thread_id][reads[i]->getIdOfCompRevRead()] = true; // comprev read is a proper suffix of another read, we remove it as well.
                }
            }
        }
    };

    vector<std::future<void> > futures(Params::THREADS - 1);
    int W = (int) ceil((double) (N - 1) / Params::THREADS);
    for (int i = 1; i < Params::THREADS; i++) {
        int a = i * W;
        int b = min((i + 1) * W - 1, N - 2);
        futures[i - 1] = std::async(std::launch::async, [&lcpFun](int a, int b, int thread_id) {
            lcpFun(a, b, thread_id);
        }, a, b, i);
    }

    lcpFun(0, min(W - 1, N - 2), 0);
    for (auto &p : futures) p.get();

    cerr << "lcp created, merging to one bitvector" << flush;
    for (int i = 0; i < markedToRemove[0].size(); i++) {
        for (int j = 1; j < Params::THREADS; j++) {
            if (markedToRemove[j][i]) markedToRemove[0][i] = true;
        }
    }
    cerr << "   merged!" << endl;

    TimeMeasurer::stopMeasurement("LCP");
    return markedToRemove[0];

}

vector<Read *> ReadPreprocess::getSortedReads() {

    TimeMeasurer::startMeasurement("SORTING");
    cerr << "Removing prefix reads - sorting" << endl;

    // THE CODE BELOW IS PARALLEL SORTING USING FIRST BUCKET SORT for first block of each biset, then usual sort for each bucket done in every thread separately.
    // creating buckets with reads
    vector<vector<Read *> > buckets(Params::THREADS);
    unsigned N = Global::READS.size();
    int MOD = (1 << 11);
    for (Read *r : Global::READS) {
        if (r == nullptr) continue;
        Bitset::TYPE bl = r->getSequence().getBlock(0);

        Bitset::TYPE bl0 = 0;
        int p = 1;
        while (p < MOD) {
            bl0 = (bl0 << 1) + (bl & 1);
            bl >>= 1;
            p <<= 1;
        }

        if (bl0 >= MOD) cerr << "bl0 >= MOD" << endl;

        bl0 &= MOD - 1; // this is just bl0 % MOD
        int ind = ((double) bl0 / MOD) * Params::THREADS;
        buckets[ind].push_back(r);
    }

    cerr << "buckets created" << endl;


    // function to sort a single bucket
    auto sortFun = [&buckets](int i) {
        sort(buckets[i].begin(), buckets[i].end(), [](Read *r1, Read *r2) {
            Bitset *b1 = &r1->getSequence();
            Bitset *b2 = &r2->getSequence();
            int m = min(b1->countBlocks(), b2->countBlocks());
            int p = 0;
            while (p < m) {
                if (b1->getBlock(p) == b2->getBlock(p)) p++;
                else {
                    int ind;
                    if (sizeof(Bitset::TYPE) == 8) ind = __builtin_ctzll(b1->getBlock(p) ^ b2->getBlock(p));
                    else ind = __builtin_ctz(b1->getBlock(p) ^ b2->getBlock(p));
                    return (*b1)[p * Bitset::BLOCK_SIZE + ind] < (*b2)[p * Bitset::BLOCK_SIZE + ind];

//                    assert( b1->getBlock(p) ^ b2->getBlock(p) != 0 ); // #TEST
                }
            }


            if (r1->size() != r2->size()) return r1->size() < r2->size();
            else return r1->getId() < r2->getId();

        });
    };



    // sorting paralelly
    vector<std::future<void> > futures(Params::THREADS - 1);
    for (int i = 1; i < Params::THREADS; i++) {
        futures[i - 1] = std::async(std::launch::async, [&sortFun](int i) {
            sortFun(i);
        }, i);
    }
    sortFun(0);
    for (auto &p : futures) p.get();

    cerr << "sorted" << endl;

    vector<Read *> reads;
    for (auto &v : buckets) {
        reads.insert(reads.end(), v.begin(), v.end());
        vector<Read *>().swap(v);
    }


    TimeMeasurer::stopMeasurement("SORTING");
//    return std::move(reads);
    return reads;
}

vector<short> ReadPreprocess::getLCP(vector<Read *> &reads) {
    cerr << "Removing prefix reads - calculating LCP" << endl;
    TimeMeasurer::startMeasurement("LCP");
    int N = reads.size();
    vector<short> lcp(N - 1, 0);

    // parallel lcp calculation - for LUX data (that is 20^6 reads) this work slower than sequential lcp.
    auto lcpFun = [&reads, &lcp](int c, int d) {
        for (int i = c; i <= d; i++) {
            int a = i;
            int b = i + 1;
            Read *r1 = reads[a];
            Read *r2 = reads[b];

            int m = min(r1->size(), r2->size());
            int l = 0;
            int BEG = 0;
            if (r1->getSequence().getBlock(0) == r2->getSequence().getBlock(0)) l = BEG = Bitset::BLOCK_SIZE / 2;

            for (int j = BEG; j < m; j++) {
                if ((*r1)[j] == (*r2)[j]) l++;
                else break;
            }
            lcp[i] = l;
        }
    };

    vector<std::future<void> > futures(Params::THREADS - 1);
    int W = (int) ceil((double) (N - 1) / Params::THREADS);
    for (int i = 1; i < Params::THREADS; i++) {
        int a = i * W;
        int b = min((i + 1) * W - 1, N - 2);
        futures[i - 1] = std::async(std::launch::async, [&lcpFun](int a, int b) {
            lcpFun(a, b);
        }, a, b);
    }

    lcpFun(0, min(W - 1, N - 2));
    for (auto &p : futures) p.get();

    TimeMeasurer::stopMeasurement("LCP");
    return lcp;
}
