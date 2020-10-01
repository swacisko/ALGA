/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: sylwester
 *
 * Created on November 23, 2018, 7:59 PM
 */

#include<iostream>
#include <experimental/filesystem>
#include <AlignmentControllers/AlignmentControllerHybrid.h>
#include <GraphSimplifiers/GraphSimplifier.h>
#include <GraphCreators/GraphCreatorKmerBased.h>
#include <GraphCreators/GraphCreatorLI.h>
#include <GraphCreators/GraphCreatorPrefSuf.h>
#include <thread>
#include <ContigCreators/ContigCreatorSinglePath.h>
#include <IO/OutputWriterNew.h>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <GraphCreators/GraphCreatorPairwiseKmerBranch.h>
#include <ContigCreators/ContigCorrector.h>
#include <assert.h>

#include "IO/InputReader.h"
#include "Global.h"
#include "IO/ReadPreprocess.h"

#include "Corrector/ReadCorrector.h"

#include "StatisticsGenerators/GenomeStatisticsCollector.h"
#include "Utils/TimeMeasurer.h"
#include "Utils/GraphVisualizer.h"

using namespace std;

void initilizeStaticData() {
    Read::priorities = VI(4);
    iota(Read::priorities.begin(), Read::priorities.end(), 0);
    Bitset::initializeStaticBlock();
}


void testConcurrency() {


    int N = 60;
    vector<Kmer> kmers(N);

    int T = 6;

    auto fun = [&kmers](int a, int b, int thr) {
        for (int i = a; i <= b; i++) {
            kmers[i].length++;
            kmers[i].hash++;
            kmers[i].indInRead *= 3;
            kmers[i].indInRead %= 1'000'001;
        }
    };

    int R = 10'000;
    for (int r = 0; r < R; r++) {
        vector<std::thread> parallelJobs;
        parallelJobs.reserve(Params::THREADS);

        int W = (int) ceil((double) N / T);
        for (int i = 1; i < T; i++) {
            int a = i * W;
            int b = min((i + 1) * W - 1, (int) N - 1);

            parallelJobs.push_back(thread([=] { fun(a, b, i); }));
        }

        fun(0, W - 1, 0);

        for (auto &p : parallelJobs) p.join();
    }

    auto kmersConcurrent = kmers;
    kmers = vector<Kmer>(N);
    for (int r = 0; r < R; r++) fun(0, N - 1, 0);

    for (int i = 0; i < N; i++) {
        assert(kmers[i].length == kmersConcurrent[i].length);
        assert(kmers[i].hash == kmersConcurrent[i].hash);
        assert(kmers[i].indInRead == kmersConcurrent[i].indInRead);
    }

    cerr << "Concurrency passed" << endl;
    exit(1);
}


class BitsetTest {
public:
    unsigned *V;
//    vector<unsigned> V;
//    deque<int> V;
    int lem;
    unsigned N;
//    char lem2;
};

class MyVec {
    char c;
    short t;
    char d;
};


int main(int argc, char **argv) {
    initilizeStaticData();

//    testConcurrency();

//    Bitset::test();
    DEBUG(sizeof(BitsetTest));
    DEBUG(sizeof(std::mutex));
    DEBUG(sizeof(MILPII));
    DEBUG(sizeof(LPII));
    DEBUG(sizeof(PII));
    DEBUG(sizeof(vector<int>));
    DEBUG(sizeof(PII[3]));
    DEBUG(sizeof(pair<unsigned, char>));
    DEBUG(sizeof(Kmer));
    DEBUG(sizeof(MyVec));
    DEBUG(sizeof(pair<short, char>));
    DEBUG(sizeof(Read));
    DEBUG(sizeof(Bitset));

//    exit(1);


    ios_base::sync_with_stdio(0);
    cin.tie(NULL);
    cout << fixed;
    cout.precision(3);
    cerr << fixed;
    cerr.precision(3);

    Params::initializeParams(argc, argv);

    const bool debug_small_test = false;
    if (debug_small_test) {
        DEBUG(std::experimental::filesystem::current_path());
        Params::THREADS = 1;
        Params::INPUT_FILE_TYPE = Params::FASTA;
        Params::fileExtension = "fasta";
        Params::TEST_NAME = "ALGA_read_test";
        Params::inStreamFilePath1 = "test_reads_1.fasta";
        Params::inStreamFilePath2 = "test_reads_2.fasta";
        Params::inStream = ifstream(Params::inStreamFilePath1);
        assert(Params::inStream.is_open());
        cin.rdbuf(Params::inStream.rdbuf());
    }


    TimeMeasurer::startMeasurement(TimeMeasurer::TOTAL_TIME);


    TimeMeasurer::startMeasurement(TimeMeasurer::INPUT_READER);
    {
        InputReader reader;
        reader.readInput();
    }
    TimeMeasurer::stopMeasurement(TimeMeasurer::INPUT_READER);
    cerr << "input read" << endl;
    MyUtils::process_mem_usage();

    vector<Read *> *READS = &Global::READS;


    int LEN = Global::calculateAvgReadLength() + Params::READ_END_TRIM_LEFT + Params::READ_END_TRIM_RIGHT;
    Params::CONTIG_MIN_OUTPUT_LENGTH = max(Params::CONTIG_MIN_OUTPUT_LENGTH, (int) (1.75 * LEN));
//    Params::CONTIG_MIN_OUTPUT_LENGTH = max( Params::CONTIG_MIN_OUTPUT_LENGTH, (int)( 1.5 * LEN) );
    Params::MAX_OFFSET_PARALLEL_PATHS = max(Params::MAX_OFFSET_PARALLEL_PATHS, (int) (1.75 * LEN));
    Params::MAX_OFFSET_DANGLING_BRANCHES = max(Params::MAX_OFFSET_DANGLING_BRANCHES, (int) (1.75 * LEN));


    if (Params::MOST_FREQUENTLY_USED_PARAMETER == -1) { // if no parameters were specified, I use default values
//        int LEN = Global::calculateAvgReadLength() + Params::READ_END_TRIM_LEFT + Params::READ_END_TRIM_RIGHT;
        int L = LEN * Params::SCALE;
//        DEBUG(Params::SCALE);
//        int RSOEMO = LEN * ( 2*Params::SCALE + 3 ) / 5;
//        int RSOEMO = LEN * ( 3*Params::SCALE + 2 ) / 5;
        int RSOEMO = LEN * (Params::SCALE + 1) / 2;

        Params::LI_KMER_LENGTH = min(2 * L / 3, 60);
        Params::KMER_LENGTH_BUCKET = min(2 * L / 3, 60);

        Params::MIN_OVERLAP_PREF_SUF = L;
        if (Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP == -1)
            Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP = RSOEMO;
        Params::MOST_FREQUENTLY_USED_PARAMETER = L;
        Params::MIN_OVERLAP_AREA = L;

    } else if (Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP == -1) {
//        int RSOEMO =  ( 50 + Params::MOST_FREQUENTLY_USED_PARAMETER/2 ) * LEN / 100;
        int RSOEMO = (Params::MOST_FREQUENTLY_USED_PARAMETER + LEN) / 2;
        Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP = RSOEMO;
    }


    Params::createOutputFileName();
    Params::writeParams();


    if (Params::CORRECT_READS) {
        ReadCorrector rc((*READS), 5, 30);
        rc.correct();
        Global::generateFasta(Params::TEST_NAME.substr(0, 8) + "_algacorrect");
        if (Params::CORRECT_READS == 2) return 0;
    }


    GenomeStatisticsCollector::addData("number of reads: ", Global::READS.size());


    if (Params::REMOVE_PREF_READS_TYPE != Params::PREF_READS_NONE) {
        ReadPreprocess prepr;
        VB prefReads = prepr.getPrefixReads();
        int cnt = 0;
        for (int i = 0; i < READS->size(); i++) {
            if (prefReads[i]) {
                Global::removeRead(i);
                cnt++;
            }
        }
        cerr << "There are " << cnt << " reads that are "
             << ((Params::REMOVE_PREF_READS_TYPE == Params::PREF_READS_ONLY_DUPLICATES) ? "duplicate"
                                                                                        : "prefix_or_suffix")
             << " of other reads - they were removed from graph creation process" << endl;
    }

    MyUtils::process_mem_usage();

    bool write_and_read = false;
    if (Params::ADD_PAIRED_READS == false) write_and_read = false;

    {
        unsigned id = 0;

        if (write_and_read) cerr << "WRITING REMAINING READS TO A SINGLE FASTA FILE" << endl;
        else cerr << "Remapping reads to the set of valid ids" << endl;
        ofstream str;
        if (write_and_read) str.open(Params::TEST_NAME + "remaining_reads.fasta");
        Global::pairedReadOffset.clear();

        unsigned back_index = 0;

        auto swap_reads = [&back_index](int i) {
            Global::READS[back_index] = Global::READS[i];
            Global::READS[back_index + 1] = Global::READS[i + 1];
            back_index += 2;
        };


        for (unsigned i = 0; i < Global::READS.size(); i += 2) {

            if (Global::READS[i] != nullptr) {
                ++id;

                assert(Global::READS[i + 1] !=
                       nullptr); // this should not happen - either both read r and its reverse complimentary are present, or neither.

                if (write_and_read) str << "read_" << id << "\n";
                if (write_and_read) str << Global::READS[i]->getSequenceAsString() << "\n";


                if ((i % 4) == 0) {

                    if (i + 2 < Global::READS.size() && Global::READS[i + 2] != nullptr) {

                        Global::pairedReadOffset.push_back(1); // for i-th read
                        Global::pairedReadOffset.push_back(1); // for reverse complimentary read

                        if (!write_and_read) swap_reads(i);

                        Global::pairedReadOffset.push_back(2); // for paired reads
                        Global::pairedReadOffset.push_back(2); // for reverse complimentary paired read

                        if (!write_and_read) swap_reads(i + 2);

                    } else {
                        Global::pairedReadOffset.push_back(0);
                        Global::pairedReadOffset.push_back(0); // for reverse complimentary read

                        if (!write_and_read) swap_reads(i);
                    }
                } else if (Global::READS[i - 2] == nullptr) {
                    Global::pairedReadOffset.push_back(0); // for this read - its paired read was removed
                    Global::pairedReadOffset.push_back(0); // for reverse complimentary read

                    if (!write_and_read) swap_reads(i);
                }
            }

            if (write_and_read && (i % 4) == 2) {
                Global::removeRead(i - 2);
                Global::removeRead(i);
            }

            if (!write_and_read && (i % 4) == 2 && 2 * id != back_index) { // just for test
                DEBUG(i);
                DEBUG(id);
                DEBUG(2 * id);
                DEBUG(back_index);
                assert(2 * id == back_index);
            }
        }
        if (write_and_read) {
            cerr << "The are " << id << " reads written to .fasta file (no complimentary reverse reads here)" << endl;
            str.flush();
            str.close();
        } else {
            DEBUG(id);
            DEBUG(2 * id);
            DEBUG(back_index);
        }

        if (write_and_read) vector<Read *>().swap(Global::READS);
        else {
            DEBUG(Global::READS.size());
            Global::READS.resize(back_index);
            vector<Read *>(Global::READS.begin(), Global::READS.end()).swap(Global::READS);
            DEBUG(Global::READS.size());
            for (int i = 0; i < Global::READS.size(); i++) Global::READS[i]->setId(i);
        }
    }

    if (write_and_read) {
        Params::inStreamFilePath1 = Params::TEST_NAME + "remaining_reads.fasta";
        Params::INPUT_FILE_TYPE = Params::FASTA;

        auto tl = Params::READ_END_TRIM_LEFT, tr = Params::READ_END_TRIM_RIGHT;
        Params::READ_END_TRIM_LEFT = Params::READ_END_TRIM_RIGHT = 0;

        Params::ADD_PAIRED_READS = 0;
//        DEBUG(Params::ADD_COMP_REV_READS);
//        Params::ADD_COMP_REV_READS = 1;
        InputReader reader;
        reader.readInput();
        Params::ADD_PAIRED_READS = 1;

        Params::READ_END_TRIM_LEFT = tl;
        Params::READ_END_TRIM_RIGHT = tr;

        DEBUG(Global::pairedReadOffset.size());
        DEBUG(Global::READS.size());
        assert(Global::READS.size() == Global::pairedReadOffset.size());
    }

    READS = &Global::READS;

    assert(Global::READS.size() == Global::pairedReadOffset.size());

    Global::GRAPH = Graph(Global::READS.size());
    Graph *G = &Global::GRAPH;

    if (Params::DESERIALIZE_GRAPH && G->deserializeGraph(Params::TEST_NAME + "_beforeSimplifier.graph"));
    else {
        TimeMeasurer::startMeasurement("GraphCreator PrefSuf");

        GraphCreator *graphCreator;
        cerr << "Creating GraphCreator" << endl;
        if (Params::ALGORITHM_IN_USE == Params::PREF_SUF_GRAPH_CREATION)
            graphCreator = new GraphCreatorPrefSuf(READS, G);
        else if (Params::USE_LI) graphCreator = new GraphCreatorLI(&Global::READS, &Global::GRAPH);


        for (int i = 0; i < READS->size(); i++) {
            if ((*READS)[i] != nullptr && (*READS)[i]->size() < Params::LI_KMER_INTERVALS + Params::LI_KMER_LENGTH) {
                graphCreator->setAlignFrom(i, false);
                graphCreator->setAlignTo(i, false);
            }
        }


        for (int i = 0; i < READS->size(); i++) {
            if ((*READS)[i] != nullptr && graphCreator->getAlignFrom(i) == false &&
                graphCreator->getAlignTo(i) == false) {
                Global::removeRead(i);
            }
        }


        int removedReads = 0;
        LL removedReadsControlSum = 0;
        for (int i = 0; i < G->size(); i++) {
            if ((*READS)[i] == nullptr) {
                removedReadsControlSum += (LL) i;
                graphCreator->setAlignFrom(i, false);
                graphCreator->setAlignTo(i, false);
                removedReads++;
            }
        }
        cerr << "Removed dispensible reads - " << removedReads << " reads were removed" << endl;
        DEBUG(removedReadsControlSum);

        graphCreator->startAlignmentGraphCreation();
        MyUtils::process_mem_usage();


        graphCreator->clear();
        delete graphCreator;
        graphCreator = 0;

        cerr << "retainingOnlySmallestOffset" << endl;
        G->retainOnlySmallestOffset();

        if (Params::SERIALIZE_GRAPH_BEFORE_SIMPLIFIER) G->serializeGraph(Params::TEST_NAME + "_beforeSimplifier.graph");

        TimeMeasurer::stopMeasurement("GraphCreator PrefSuf");
    }

    MyUtils::process_mem_usage();

//    bool usePkbSupplement = ( Params::TPN.find( "suppl" ) != string::npos );
    bool usePkbSupplement = Params::USE_GRAPH_CREATOR_SUPPLEMENT;
//    bool usePkbSupplement = true;
    if (usePkbSupplement) {
        TimeMeasurer::startMeasurement("GraphCreator PKB Supplement");

        cerr << "Before supplement, G has " << G->countEdges() << " edges" << endl;
        GraphCreator *graphCreator = new GraphCreatorLI(READS, G);

        VI *inDeg = G->getInDegrees();

        for (int i = 0; i < G->size(); i++) {
            graphCreator->setAlignFrom(i, false);
            graphCreator->setAlignTo(i, false);

            if ((*inDeg)[i] == 0 && (*G)[i].size() > 0) graphCreator->setAlignTo(i, true);
            if ((*inDeg)[i] > 0 && (*G)[i].size() == 0) graphCreator->setAlignFrom(i, true);

            if (Params::USE_GRAPH_CREATOR_SUPPLEMENT == 1 && (*inDeg)[i] > 0 && (*G)[i].size() > 0 &&
                ((*inDeg)[i] >= 2 || (*G)[i].size() >= 2)) {
                graphCreator->setAlignFrom(i, true);
                graphCreator->setAlignTo(i, true);
            }
        }

        Params::MIN_OVERLAP_AREA = Global::calculateAvgReadLength() -
                                   (Params::USE_GRAPH_CREATOR_SUPPLEMENT ? 10 : 0); // the latter is for GCPS

        Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT = (Params::USE_GRAPH_CREATOR_SUPPLEMENT ? 10
                                                                                            : 1);    // the latter is for GCPS
        Params::MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = (Params::USE_GRAPH_CREATOR_SUPPLEMENT ? 95
                                                                                          : 99);   // the latter is for GCPS

        Params::LI_KMER_INTERVALS = 6;
        Params::LI_KMER_LENGTH = 35;
        graphCreator->startAlignmentGraphCreation();


        cerr << "retainingOnlySmallestOffset" << endl;
        G->retainOnlySmallestOffset();

        delete inDeg;
        inDeg = 0;

        delete graphCreator;
        graphCreator = 0;
        cerr << "After supplement G has " << G->countEdges() << " edges" << endl;

        TimeMeasurer::stopMeasurement("GraphCreator PKB Supplement");
    }


    G->pruneGraph();
    Global::removeIsolatedReads();

    {
        cerr << "Before first simplifier graph has " << G->countEdges() << " edges" << endl;

        GraphSimplifier simplifier(Global::GRAPH, *READS);


        if (Params::DESERIALIZE_GRAPH &&
            G->deserializeGraph(Params::TEST_NAME + "_mopp" + to_string(Params::MAX_OFFSET_PARALLEL_PATHS)
                                + "_modb" + to_string(Params::MAX_OFFSET_DANGLING_BRANCHES) + "_rsoe" +
                                to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP) + "-" +
                                to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN) +
                                "_afterSimplifier.graph")) {

        } else {
            simplifier.simplifyGraphOld();

            if (Params::SERIALIZE_GRAPH_AFTER_SIMPLIFIER)
                G->serializeGraph(Params::TEST_NAME + "_mopp" + to_string(Params::MAX_OFFSET_PARALLEL_PATHS)
                                  + "_modb" + to_string(Params::MAX_OFFSET_DANGLING_BRANCHES) + "_rsoe" +
                                  to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP) + "-" +
                                  to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN) +
                                  "_afterSimplifier.graph");
            cerr << "Memory usage after simplifyGraphOld, before pruning graph and raeds" << endl;
            MyUtils::process_mem_usage();
        }

        G->pruneGraph();
//        cerr << "Removing isolated reads" << endl;
        Global::removeIsolatedReads();

//        cerr << "Creating contracted edges vector" << endl;
        G->createContractedEdgesVector();


        cerr << endl << "AFTER SIMPLIFIER, before contracting paths" << endl;
        MyUtils::process_mem_usage();


        for (int x = 0; x < 2; x++) {
            G->retainOnlySmallestOffset();
            simplifier.simplifyGraph();
        }

        cerr << endl << "AFTER SIMPLIFIER - CONTRACTING PATHS" << endl;
        MyUtils::process_mem_usage();

    } // end of simplifier


    G->retainOnlySmallestOffset();


    vector<Contig *> contigs;
    {
        ContigCreator *ctgCreator = new ContigCreatorSinglePath(G, *READS);
        contigs = ctgCreator->getAllContigs();
        delete ctgCreator;
        ctgCreator = 0;

        cerr << "Contigs created" << endl;
        MyUtils::process_mem_usage();
    }




//    {
//        GraphVisualizer gviz;
//        gviz.writeWholeGraph(G, *READS, Params::TEST_NAME + "_graph_gviz");
//        gviz.writeInGraphvizFormat(G, contigs);
//    }


    cerr << "Before filtering, there are " << contigs.size() << " contigs" << endl;
    MyUtils::process_mem_usage();
    contigs = OutputWriterNew(G, contigs).filterContigs();
    cerr << endl << "After filtering, there are " << contigs.size() << " contigs" << endl;
    MyUtils::process_mem_usage();

    cerr << "Contigs length check:" << endl;
    for (auto t : contigs) {
        if (t->size() != t->getSequenceAsString().size()) {
            DEBUG(t->size());
            DEBUG(t->getSequenceAsString().size());
            exit(1);
        }
    }


    for (auto ctg : contigs) {
        for (auto pair : ctg->getContainedReads()) {
            if (pair.first == nullptr) {
                cerr << "pair.first == nullptr" << endl;
                exit(1);
            } else if (pair.first->getId() < 0 || pair.first->getId() >= G->size()) {
                cerr << "pair.first->getId() = " << pair.first->getId() << endl;
            }
        }
    }



//        bool extendContigs = true;
    bool extendContigs = Params::TPN.find("extend") != string::npos;
    if (extendContigs) {

        vector<Read *> newReads;
        for (auto t : contigs) newReads.push_back(t);

        int cnt = 0;
        for (auto t : contigs)
            newReads.push_back(
                    new Read(cnt++, MyUtils::getComplimentaryString(MyUtils::getReverse(t->getSequenceAsString()))));

        cnt = 0;
        for (auto t : newReads) t->setId(cnt++);

        Graph *newGraph = new Graph(newReads.size());
        GraphCreatorPrefSuf gcps(&newReads, newGraph);

        int THRESHOLD = 30;
        Params::MIN_OVERLAP_PREF_SUF = THRESHOLD;
        Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP = THRESHOLD;

        gcps.startAlignmentGraphCreation();


        int M = newReads.size() / 2;

//            GraphVisualizer gviz;
//            gviz.writeWholeGraph(newGraph,newReads, "ALGA_not_extended");

        auto countPEConnections = [&contigs, &M](int a, int b) {
            vector<pair<Read *, int>> ce1, ce2;

            if (a >= M) a -= M;
            ce1 = contigs[a]->getContainedReads();

            if (b >= M) b -= M;
            ce2 = contigs[b]->getContainedReads();

            unordered_set<int> zb1;

            int DST = 500;
            for (int i = 0; i < min(DST, (int) ce1.size()); i++) {
                zb1.insert(ce1[i].first->getId());
                zb1.insert(ce1[ce1.size() - 1 - i].first->getId());
            }

            int cnt = 0;
            for (int i = 0; i < min(DST, (int) ce2.size()); i++) {
                if (zb1.count(Read::getIdOfPairedRead(ce2[i].first->getIdOfCompRevRead()))
                    || zb1.count(Read::getIdOfPairedRead(ce2[ce2.size() - 1 - i].first->getIdOfCompRevRead()))) {
                    cnt++;
                }
            }

            return cnt;
        };


        int PE_THRESHOLD = 5;

        VPII toRemove;
        for (int i = 0; i < 2 * M; i++) {
            for (PII neigh : (*newGraph)[i]) {
                int d = neigh.first;
                int cnt = countPEConnections(i, d);
                if (cnt < PE_THRESHOLD) {
                    toRemove.push_back({i, d});
//                        cerr << "edge (" << i << "," << d << ") has only " << cnt << " connections - removing!" << endl;
                } else {
//                        cerr << endl;
//                        cerr << "\t\tBINGO!" << "Edge (" << i << "," << d << ") has " << cnt << " connections!" << endl;
//                        cerr << endl;
                }
            }
        }

        for (PII e : toRemove) newGraph->removeDirectedEdge(e.first, e.second);


//            gviz.writeWholeGraph(newGraph,newReads, "ALGA_removed_PE");


        Graph newGraphRev = newGraph->getReverseGraph();
        unordered_set<int> contigsToRemove;

        auto contractDisjointPaths = [&M, &newGraph, &newGraphRev, &newReads, &contigsToRemove]() {
            VVI paths;

            for (int i = 0; i < 2 * M; i++) {

//                    cerr << "i = " << i << endl;
                if ((*newGraph)[i].size() != 1 || newGraphRev[i].size() > 1) continue;

//                    cerr << "\tconsidering i = " << i << endl;

                for (auto neigh : (*newGraph)[i]) {
                    int d = neigh.first;

//                        cerr << "\td = " << d << endl;

                    if ((*newGraph)[d].size() > 1 || newGraphRev[d].size() > 1) continue;

                    // now i extend contig i so that it contains also contig d

                    int offset = neigh.second;
                    int overlap = Read::calculateReadOverlap(newReads[i], newReads[d], offset);

                    int x = -1;

                    bool hasSuccessor = (*newGraph)[d].size() > 0;
                    if (hasSuccessor) x = (*newGraph)[d][0].first;
                    int offx;
                    if (hasSuccessor) offx = (*newGraph)[d][0].second;


//                        cerr << "i = " << i << "   d = " << d << "    x = " << x << endl;
//                        cerr << "contig[i].size() = " << newReads[i]->size() << "   contigs[d].size() = " << newReads[d]->size() << endl;
//                        cerr << "offset(i,d) = " << offset << endl;
//                        if(hasSuccessor)cerr << "offset(d,x) = " << offx << endl << endl;

                    string s = newReads[i]->getSequenceAsString();
                    s = s.substr(0, s.size() - overlap);
                    s += newReads[d]->getSequenceAsString();
                    newReads[i]->modifySequence(s);


                    contigsToRemove.insert(d);


                    if (hasSuccessor) {
                        (*newGraph).pushDirectedEdge(i, x, offset + offx);
                        newGraphRev.pushDirectedEdge(x, i, offset + offx);

                        (*newGraph).removeDirectedEdge(d, x);
                        newGraphRev.removeDirectedEdge(x, d);
                    }

                    (*newGraph).removeDirectedEdge(i, d);
                    newGraphRev.removeDirectedEdge(d, i);

                }
            }
        };


        for (int i = 0; i < 3; i++) contractDisjointPaths();

        int p = M - 1;
        for (int i = p; i >= 0; i--) {
            if (contigsToRemove.count(i)) {
                swap(contigs[i], contigs.back());
                contigs.pop_back();
                (*newGraph)[i].clear();
                for (auto p : newGraphRev[i]) (*newGraph).removeDirectedEdge(p.first, i);
            }
        }

//            gviz.writeWholeGraph(newGraph,newReads, "ALGA_extended");
    }


    bool trimContigs = true;
//        bool trimContigs = Params::TPN.find( "trim" ) != string::npos;
    if (trimContigs) {


        vector<Read *> newReads;
        for (auto t : contigs) newReads.push_back(t);

        int cnt = 0;
        for (auto t : contigs)
            newReads.push_back(
                    new Read(cnt++, MyUtils::getComplimentaryString(MyUtils::getReverse(t->getSequenceAsString()))));

        cnt = 0;
        for (auto t : newReads) t->setId(cnt++);

        Graph *newGraph = new Graph(newReads.size());
        GraphCreatorPrefSuf gcps(&newReads, newGraph);

        int THRESHOLD = 25;
        Params::MIN_OVERLAP_PREF_SUF = THRESHOLD;
        Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP = THRESHOLD;

        gcps.startAlignmentGraphCreation();


        for (auto &ctg : contigs) {
            for (auto &pair : ctg->getContainedReads()) {
                if (pair.first == nullptr) {
                    cerr << "pair.first == nullptr" << endl;
                    exit(1);
                } else if (pair.first->getId() < 0 || pair.first->getId() >= G->size()) {
                    cerr << "pair.first->getId() = " << pair.first->getId() << endl;
                    exit(1);
                }
            }
        }


//            GraphVisualizer gviz;
//            gviz.writeWholeGraph(newGraph,newReads);



        cerr << "Trimming!" << endl;

        VI trimLeft(newReads.size() / 2, 0);
        VI trimRight(newReads.size() / 2, 0);

        int M = newReads.size() / 2;
        int avg_read_length = Global::calculateAvgReadLength();

        for (int i = 0; i < newReads.size(); i++) {

            for (PII neigh : (*newGraph)[i]) {
                int d = neigh.first;
                int offset = neigh.second;
                int overlap = newReads[i]->size() - offset;

                if (overlap >= 2 * avg_read_length) {
                    cerr << "Found pair with overlap " << overlap << endl;
                }

                if (i < M && d < M) trimLeft[d] = max(trimLeft[d], overlap);
//                    else{
//                        if( i < M && newReads[i]->size() < newReads[d]->size() ) trimRight[i] = max( trimRight[i], overlap );
//                        if( d < M && newReads[d]->size() < newReads[i]->size() ) trimLeft[d] = max( trimLeft[d], overlap );
//                    }

            }
        }


        for (int i = 0; i < M; i++) {


            Contig *t = (Contig *) newReads[i];
            string s = t->getSequenceAsString();
            if (trimLeft[i] + trimRight[i] + 10 < (int) s.size())
                s = s.substr(trimLeft[i], max(1, (int) s.size() - trimLeft[i] - trimRight[i]));
            else s = "CCCC";

            Read *r = new Read(0, s);
            t->modifySequence(s);
            delete r;

        }

        for (int i = M; i < 2 * M; i++) {
            delete newReads[i];
            newReads[i] = nullptr;
        }

        delete newGraph;
        newGraph = nullptr;

//            cerr << "Trimmed" << endl;


        Params::ADD_PAIRED_READS = Params::ADD_COMP_REV_READS = 0;
        cerr << "Leaving trimming section" << endl;
    }


    TimeMeasurer::startMeasurement(TimeMeasurer::OUTPUT_WRITER);

//        DEBUG(contigs.size());
    OutputWriterNew writer(G, contigs);
//        writer.writeContigs();
    writer.writeContigsNoFilter(contigs);
    TimeMeasurer::stopMeasurement(TimeMeasurer::OUTPUT_WRITER);


    TimeMeasurer::stopMeasurement(TimeMeasurer::TOTAL_TIME);

    //************************  WRITING STATISTICS AND MEASUREMENTS

    Params::WRITE_STATISTICS = 1;
    Params::writeParams();
    TimeMeasurer::writeAllMeasurements();
    GenomeStatisticsCollector::writeTestStatistics();

    cerr << endl << "READS WITH N STATISTICS:" << endl;
    // for(auto a : reader.getReadsWithNSizes()) cerr << a << " " << flush;
    //    StatisticsGenerator::writeAllStatistics( reader.getReadsWithNSizes() );

    StatisticsGeneratorBigData::writeAllStatistics();

    Params::WRITE_STATISTICS = 1;
    cerr << endl << "CONTIG-LENGTHS STATISTICS:" << endl;
    StatisticsGenerator::writeAllStatistics(writer.getContigsLengths());


    // clearing dynamically allocated memory
    {
        for (Contig *c : contigs) {
            if (c != nullptr) {
                delete c;
                c = nullptr;
            }
        }

        for (int i = 0; i < Global::READS.size(); i++) {
            if (Global::READS[i] != nullptr) Global::removeRead(i);
        }


    }

    return 0;

}

