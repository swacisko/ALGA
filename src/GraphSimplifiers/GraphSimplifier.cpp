/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   GraphSimplifier.cpp
 * Author: sylwester
 *
 * Created on November 25, 2018, 12:42 AM
 */

#include <AlignmentControllers/AlignmentControllerLCS.h>
#include <AlignmentControllers/AlignmentControllerHybrid.h>
#include <GraphSimplifiers/GraphSimplifier.h>
#include <Utils/TimeMeasurer.h>
#include <StatisticsGenerators/GenomeStatisticsCollector.h>
#include <AlignmentControllers/AlignmentControllerLowErrorRate.h>
#include <thread>
#include <queue>
#include <Utils/WorkloadManager.h>

#include "GraphSimplifiers/GraphSimplifier.h"
#include "Global.h"
#include "DataStructures/FAU.h"
#include "GraphCreators/GraphCreator.h"
#include "GraphCreators/GraphCreatorPrefSuf.h"

GraphSimplifier::GraphSimplifier(Graph &G, vector<Read *> &reads) {
    this->G = &G;
    this->reads = &reads;
    was = VVB(Params::THREADS, VB(this->G->size(), false));
    edgesToRemove = VVPII(Params::THREADS);
}

GraphSimplifier::GraphSimplifier(const GraphSimplifier &orig) {
}

GraphSimplifier::~GraphSimplifier() {
    clear();
}

void GraphSimplifier::clear() {
    VVB().swap(was);
    VVPII().swap(edgesToRemove);
}

void GraphSimplifier::simplifyGraph() {
    TimeMeasurer::startMeasurement(TimeMeasurer::GRAPH_SIMPLIFIER);
    cerr << endl
         << "**************************************************************************** BEGINNING OF GRAPH SIMPLIFIER"
         << endl;


    int REMOVE_ITERATIONS = 1;

    for (int i = 0; i < REMOVE_ITERATIONS; i++) {
        cerr << "******************************************* ITERATION " << i + 1 << endl;

        cerr << endl << "****** BEFORE CUTTING METRIC TRIANGLES" << endl;


        cutNonAndWeaklyMetricTriangles();

        cerr << endl << "****** AFTER CUTTING METRIC TRIANGLES, BEFORE CONTRACTING PATHS" << endl;

        bool contractionDone = contractPathNodes();

        // this here is to iterate as long as any contraction is done
        if (i == REMOVE_ITERATIONS - 1 && contractionDone) REMOVE_ITERATIONS++;
    }


    TimeMeasurer::stopMeasurement(TimeMeasurer::GRAPH_SIMPLIFIER);


    cerr << endl
         << "**************************************************************************** END OF GRAPH SIMPLIFIER"
         << endl;

}


void GraphSimplifier::simplifyGraphOld() {
    TimeMeasurer::startMeasurement(TimeMeasurer::GRAPH_SIMPLIFIER);
    cerr << endl
         << "**************************************************************************** BEGINNING OF GRAPH SIMPLIFIER"
         << endl;


    int min_overlap = Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP;
    int number_to_retain = Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN;

    G->sortEdgesByIncreasingOffset();
    if (Params::ALGORITHM_IN_USE != Params::PREF_SUF_GRAPH_CREATION) {
        removeSmallOverlapEdges(min_overlap, number_to_retain);
        G->pruneGraph();
    }
//    removeSmallOverlapEdges(min_overlap, number_to_retain); G->pruneGraph();  // this should not make any change in case of GCPS graph creation

//    Global::checkOLCGraphCorrectness(G,reads);


    G->sortEdgesByIncreasingOffset();


    if (Params::ALGORITHM_IN_USE != Params::PREF_SUF_GRAPH_CREATION ||
        Params::REMOVE_PREF_READS_TYPE != Params::PREF_READS_ALL_PREFIX_READS)
        while (mergeLength0Edges()) {} // merge until impossible to merge more
//    while (mergeLength0Edges()) {} // this should not make any change in case of GCPS graph creation


    G->sortEdgesByIncreasingOffset();

    cutNonAndWeaklyMetricTriangles();
    Global::removeIsolatedReads();


    for (int i = 0; i < Params::THREADS; i++) VPII().swap(edgesToRemove[i]);


    G->pruneGraph();

    int REMOVE_ITERATIONS = 1;


    Global::calculateAvgReadLength();


    removeShortParallelPaths((int) ((double) (Params::MAX_OFFSET_PARALLEL_PATHS * Global::AVG_READ_LENGTH) /
                                    (float) 100)); // this should remove all the triangles as well. All parallel paths with length no more than 1/4 readlength will be removed


    G->pruneGraph();


    Global::removeIsolatedReads();

    G->retainOnlySmallestOffset();


    for (int i = 0; i < REMOVE_ITERATIONS; i++) {
        cerr << endl << endl
             << "******************************************************************** GraphSimplifier REMOVE_ITERATION = "
             << (i + 1) << endl;

        cerr << "\rGraphSimplifier removing dangling branches             " << flush;

        int removed = 0;
        removed += removeDanglingBranches(
                (int) (Params::MAX_OFFSET_DANGLING_BRANCHES * Global::AVG_READ_LENGTH / (float) 100));

        cerr << "\rGraphSimplifier removing upper dangling branches             " << flush;
        removed += removeDanglingUpperBranches(
                (int) (Params::MAX_OFFSET_DANGLING_BRANCHES * Global::AVG_READ_LENGTH / (float) 100));


        if (i == REMOVE_ITERATIONS - 1 && removed > 0) REMOVE_ITERATIONS++;

        if (i >= 15 && removed <= 30)
            break; // just to prevent very long loops that sometimes occur in the graph and make absolutely no influence on the quality of the results.

    }

    Global::removeIsolatedReads();

    TimeMeasurer::stopMeasurement(TimeMeasurer::GRAPH_SIMPLIFIER);


    cerr << endl
         << "**************************************************************************** END OF GRAPH SIMPLIFIER"
         << endl;

}

void GraphSimplifier::cutNonAndWeaklyMetricTriangles() {

    TimeMeasurer::startMeasurement("GraphSimplifier_cutNonAndWeaklyMetricTriangles");

    for (auto &p : edgesToRemove) p.clear();
    vector<std::thread> parallelJobs;
    parallelJobs.reserve(Params::THREADS);

    int W = (int) ceil((double) G->size() / Params::THREADS);
    for (int i = 1; i < Params::THREADS; i++) {
        int a = min(i * W, G->size() - 1);
        int b = min((i + 1) * W - 1, G->size() - 1);

        parallelJobs.push_back(thread([=] { removeNonAndWeaklyMetricTrianglesJobAddToRemove(a, b, i); }));
    }

    removeNonAndWeaklyMetricTrianglesJobAddToRemove(0, W - 1, 0);

    for (auto &p : parallelJobs) p.join();


    int totalEdgesToRemove = 0;
    for (auto &p : edgesToRemove) totalEdgesToRemove += p.size();
    cerr << endl << "There are " << totalEdgesToRemove << " edges to remove in cutNonAndWeaklyMetricTriangles" << endl;

    parallelJobs.clear();
    for (int i = 1; i < Params::THREADS; i++) {
        int a = min(i * W, G->size() - 1);
        int b = min((i + 1) * W - 1, G->size() - 1);

        parallelJobs.push_back(thread([=] { removeNonAndWeaklyMetricTrianglesJobRemoveEdges(a, b, i); }));
    }

    removeNonAndWeaklyMetricTrianglesJobRemoveEdges(0, W - 1, 0);

    for (auto &p : parallelJobs) p.join();


    TimeMeasurer::stopMeasurement("GraphSimplifier_cutNonAndWeaklyMetricTriangles");

}

void GraphSimplifier::removeNonAndWeaklyMetricTrianglesJobRemoveEdges(int a, int b, int thread_id) {
    int processed = 0;
    int progressCounter = 0;
    for (auto a : edgesToRemove[thread_id]) {
        G->removeDirectedEdge(a.first, a.second);
        if (thread_id == 0)
            MyUtils::writeProgress(++processed, edgesToRemove[thread_id].size(), progressCounter,
                                   "cutting triangles progress", 1);
    }

    VPII().swap(edgesToRemove[thread_id]);
}


void GraphSimplifier::removeNonAndWeaklyMetricTrianglesJobAddToRemove(int a, int b, int thread_id) {

    int trulyMetricTriangles = 0;
    int nonAndWeaklyMetricTriangles = 0;
    int progressCounter = 0;

    unordered_map<int, int> dst;

    for (int i = a; i <= b; i++) {

        for (int k = 0; k < (*G)[i].size(); k++) {
            int a = (*G)[i][k].first;

            for (int j = 0; j < (*G)[a].size(); j++) {
                int b = (*G)[a][j].first;

                if (dst.count(b) == 0) dst[b] = G->getWeight(i, k) + G->getWeight(a, j);
                else dst[b] = min(dst[b], G->getWeight(i, k) + G->getWeight(a, j));
            }
        }


        for (int k = 0; k < (*G)[i].size(); k++) {

            int b = (*G)[i][k].first;

            if (G->getWeight(i, k) > Params::MAX_OFFSET_PARALLEL_PATHS) {
                continue; /* i do not remove long paths, even if they are in nonmetric triangle*/
            }

//            if( dst[thread_id][b] != -1 && ( dst[thread_id][b] <= G->getWeight(i,k) ) ){

//            if (dst.count(b) /*&& ( dst[b] <= G->getWeight(i,k) )*/ ) {

//            if (dst.count(b) && ( dst[b] <= G->getWeight(i,k) ) ) { // #TEST
            if (dst.count(b) && (dst[b] == G->getWeight(i, k))) { // #TEST, equals distances
//            if (dst.count(b) && ( abs( (int)dst[b] - (int)G->getWeight(i,k) ) <= 1 ) ) { // #TEST
//            if (dst.count(b) && ( dst[b] >= G->getWeight(i,k) * 0.95 && dst[b] <= G->getWeight(i,k) ) * 1.05 ) { // #TEST

                edgesToRemove[thread_id].emplace_back(i, b);
                nonAndWeaklyMetricTriangles++;
//            }else if( dst[thread_id][b] != -1 /*&& Params::REMOVE_TRULY_METRIC_TRIANGLES*/ ){
            } else if (dst.count(b) != 0 /*&& Params::REMOVE_TRULY_METRIC_TRIANGLES*/ ) {
                trulyMetricTriangles++;
            }
        }

//        for( int k=0; k< (*G)[i].size(); k++ ){ // equivalent to dst.clear()
//            int a = (*G)[i][k].first;
//            for( int j=0; j< (*G)[a].size(); j++ ){
//                int b = (*G)[a][j].first;
//                dst[thread_id][b] = -1;
//            }
//        }
        dst.clear();

        if (thread_id == 0)
            MyUtils::writeProgress(i + 1, G->size() / Params::THREADS, progressCounter,
                                   "adding triangle-edges to cut progress", 1);
    }

    if (thread_id == 0)
        cerr << "\tNonAndWeaklyMetric triangles cut in thread " << thread_id << ":  " << nonAndWeaklyMetricTriangles
             << " w eakly metric triangles cut, " << trulyMetricTriangles << " truly metric triangles present" << endl;
}


void GraphSimplifier::removeShortParallelPaths(int maxOffset) {
    TimeMeasurer::startMeasurement("GraphSimplifier_removeShortParallelPaths");
    cerr << endl << "Removing short parallel paths" << endl;

    const bool USE_WORKLOAD_MANAGER = true;
    if (!USE_WORKLOAD_MANAGER) { // original version
        vector<std::thread> parallelJobs;
        parallelJobs.reserve(Params::THREADS);

        for (int i = 1; i < Params::THREADS; i++) {
            parallelJobs.push_back(thread([=] {
                removeShortParallelPathsJob(0, G->size() - 1, maxOffset, i);
            })); // i should pass the whole scope,, thread job is done based on fau.Find value
        }

        removeShortParallelPathsJob(0, G->size() - 1, maxOffset, 0);

        for (auto &p : parallelJobs) p.join();

//        for (int i = 0; i < G->size(); i++) (*G)[i].shrink_to_fit();
    } else { // version with workload manger
        VI pathsConsidered(Params::THREADS, 0);
        VI blocksConsidered(Params::THREADS, 0);

        unsigned blocks = 50 * Params::THREADS;
        WorkloadManager::parallelBlockExecution(0, G->size() - 1, blocks, Params::THREADS,
                                                [=, &pathsConsidered, &blocksConsidered](unsigned a, unsigned b,
                                                                                         unsigned thread_id) {

                                                    for (int i = a; i <= b; i++) {
                                                        G->lockNode(i);
                                                        bool condition = ((*G)[i].size() >= 2);
                                                        G->unlockNode(i);

                                                        if (condition) {
                                                            tryToRemoveShortPathsMST(i, maxOffset, thread_id);
                                                            pathsConsidered[thread_id]++;
                                                        }
                                                    }
             blocksConsidered[thread_id]++;
                                                    if (thread_id == 0)
                                                        cerr << "\rThread 0 finished " << blocksConsidered[0]
                                                             << " block out of average "
                                                             << (double) blocks / Params::THREADS << flush;
         });

         for( int i=0; i<pathsConsidered.size(); i++){
             cerr << "There were " << pathsConsidered[i] << " paths and " << blocksConsidered[i] << " blocks considered in thread " << i << endl;
         }
    }

    cerr << "\tShortParallelPaths removed" << endl;
    TimeMeasurer::stopMeasurement("GraphSimplifier_removeShortParallelPaths");
}

void GraphSimplifier::removeShortParallelPathsJob(int a, int b, int maxOffset, int thread_id) {
    int progressCounter = 0;
    int pathsConsidered = 0;
    for (int i = a; i <= b; i++) {
        if (Params::REMOVE_SHORT_PARALLEL_PATHS_MST_FAST == true && i % Params::THREADS ==
                                                                    thread_id) { // THIS VERSION HERE IS FOR 'PERFECT' PARALLEL removeParallelPathsMST methodwith G->lockNode() in use.

            if ((*G)[i].size() >= 2 /*|| (*inDeg)[i] >= 2 */) {
                tryToRemoveShortPathsMST(i, maxOffset, thread_id);
                pathsConsidered++;
            }

            if (thread_id == 0)
                MyUtils::writeProgress(i - a + 1, b - a + 1, progressCounter,
                                       "removeShortParallelPaths progress in thread 0", 10);

        }
    }

    cerr << "There were " << pathsConsidered
         << " paths considered in removeShortParallelPaths in thread " << thread_id << endl;
}


void GraphSimplifier::tryToRemoveShortPathsMST(int beg, int maxOffset, int thread_id) {

    vector<pair<pair<int, int>, int> > edges;  // each entry is an edge of the form ( (a,b), offset )
    VI neigh(1, beg);

    unordered_map<int, int> dst, par;

    dst[beg] = 0;
    for (int i = 0; i < neigh.size(); i++) {
        int a = neigh[i];

        if (was[thread_id][a] || dst[a] > maxOffset) continue;

        was[thread_id][a] = true;

        G->lockNode(a);
        for (int k = 0; k < (*G)[a].size(); k++) {
            int b = (*G)[a][k].first;
            int offset = G->getWeight(a, k);

            if (dst.count(b) > 0 && dst[b] < dst[a] + offset) continue;


            /*if( b == beg ){ // this here is to avoid disjoining graph when we try to remove parallel paths in areas such as Short Tandem Repeats
//                cerr << "b == beg in RSPP" << endl;
//                exit(1);
//                for( auto &a : edges ) par[thread_id][ a.first.second ] = -1; // equivalent to par.clear();
                 par.clear();

                for(auto &a : edges) was[thread_id][ a.first.second ] = false; // equivalent to was.clear()

//                for( int j=0; j<neigh.size(); j++ ){ // equivalent to dst.clear();
//                    int a = neigh[j];
//                    dst[thread_id][a] = -1;
//                }
                 dst.clear();

                neigh.clear();
                G->unlockNode(a);
                return;
            }*/



            dst[b] = dst[a] + offset;

            edges.push_back({{a, b}, offset});
            neigh.push_back(b);

        }
        G->unlockNode(a);
    }

    for (auto &a : edges) {
        G->lockNode(a.first.first);
        G->removeDirectedEdge(a.first.first, a.first.second); // removing the edge from the graph   // SLOW VERSION
        G->unlockNode(a.first.first);
    }


    sort(edges.begin(), edges.end(), [](auto a, auto b) {
        if (a.second != b.second) return a.second < b.second;
        else return a.first < b.first;
    }); // sorted by MINIMAL OFFSET for directed version of MST.




    for (auto &a : neigh) was[thread_id][a] = false; // equivalent to was.clear()

    for (auto &a : edges) { // HERE I ADD EDGES ACCORDING TO MST METHOD. A NODE CANNOT BE AN END OF MORE THAN ONE EDGE (from edges) AFTER THIS.
        if (was[thread_id][a.first.second]) continue;

        G->lockNode(a.first.first);
        G->pushDirectedEdge(a.first.first, a.first.second, a.second);  // FAST VERSION
        G->unlockNode(a.first.first);

        par[a.first.second] = a.first.first;

        was[thread_id][a.first.second] = true;
    }

    par.clear();

    for (auto &a : edges) was[thread_id][a.first.second] = false; // equivalent to was.clear()

    dst.clear();

    neigh.clear();
}


int GraphSimplifier::removeDanglingBranches(int maxOffset) {
    TimeMeasurer::startMeasurement("GraphSimplifier_removeDanglingBranches");


    VVPII edgesToRemove(Params::THREADS);
    int branchesRemoved = 0;

    const bool USE_WORKLOAD_MANAGER = true;
    if (!USE_WORKLOAD_MANAGER) {

        vector<std::future<int> > futures(Params::THREADS - 1);

        int W = (int) ceil((double) G->size() / Params::THREADS);
        for (int i = 1; i < Params::THREADS; i++) {
            int a = i * W;
            int b = min((i + 1) * W - 1, (int) G->size() - 1);

            futures[i - 1] = std::async(std::launch::async, [=, &edgesToRemove](int a, int b, int i) {
                int branchesRemoved = 0;
                for (int j = a; j <= b; j++) {
                    if ((*G)[j].size() >= 2) {
                        branchesRemoved += removeDanglingBranchesFromNode(j, maxOffset, edgesToRemove[i], i);
                    }
                }

                return branchesRemoved;
            }, a, b, i);

        }


        for (int j = 0; j <= W - 1; j++) {
            if ((*G)[j].size() >= 2) {
                removeDanglingBranchesFromNode(j, maxOffset, edgesToRemove[0], 0);
            }
        }

        for (auto &p : futures) p.get();

        VPII toRemove;
        for (auto &v : edgesToRemove) {
            toRemove.insert(toRemove.end(), v.begin(), v.end());
            VPII().swap(v);
        }

        sort(toRemove.begin(), toRemove.end());
        toRemove.resize(unique(toRemove.begin(), toRemove.end()) - toRemove.begin());

        int branchesRemoved = 0;
        int progressCounter = 0;
        for (auto e : toRemove) {
            branchesRemoved++;
            if (G->removeDirectedEdge(e.first, e.second) == false) {
                cerr << "Trying to remove an edge that is not present in the graph! In removeDanglingBranches()"
                     << endl;
            }
            MyUtils::writeProgress(progressCounter + 1, G->size(), progressCounter,
                                   "RemoveDanglingBranches, removing edges", 1);
        }
        toRemove.clear();
    } else {

        VI threadBranchesRemoved(Params::THREADS, 0);

        int blocks = 10 * Params::THREADS;
        WorkloadManager::parallelBlockExecution(0, G->size() - 1, blocks, Params::THREADS,
                                                [=, &edgesToRemove](unsigned a, unsigned b, unsigned id) {
                                                    for (unsigned j = a; j <= b; j++) {
                                                        if ((*G)[j].size() >= 2) {
                                                            removeDanglingBranchesFromNode(j, maxOffset,
                                                                                           edgesToRemove[id], id);
                                                        }
                                                    }
                                                });

        { // sequential sort and parallel removal of branches
            VPII toRemove;
            for (auto &v : edgesToRemove) {
                toRemove.insert(toRemove.end(), v.begin(), v.end());
                VPII().swap(v);
            }

            sort(toRemove.begin(), toRemove.end());
            toRemove.resize(unique(toRemove.begin(), toRemove.end()) - toRemove.begin());

            /*branchesRemoved = 0;
            int progressCounter = 0;
            for (auto e : toRemove) {
                branchesRemoved++;
                if (G->removeDirectedEdge(e.first, e.second) == false) {
                    cerr << "Trying to remove an edge that is not present in the graph! In removeDanglingBranches()"
                         << endl;
                }
                MyUtils::writeProgress(progressCounter + 1, G->size(), progressCounter,
                                       "RemoveDanglingBranches, removing edges", 1);
            }*/

            // removing edges parallelly.
            if (!toRemove.empty()) {
                blocks = 3 * Params::THREADS;
                random_shuffle(toRemove.begin(),
                               toRemove.end()); // this should decrease the number of lock collisions in parallel removal of edges
                WorkloadManager::parallelBlockExecution(0, toRemove.size() - 1, blocks, Params::THREADS,
                                                        [=, &toRemove, &threadBranchesRemoved](unsigned a, unsigned b,
                                                                                               unsigned id) {
                                                            for (unsigned i = a; i <= b; i++) {
                                                                PII *e = &toRemove[i];
                                                                G->lockNode(e->first);
                                                                if (G->removeDirectedEdge(e->first,
                                                                                          e->second))
                                                                    threadBranchesRemoved[id]++;
                                                                G->unlockNode(e->first);
                                                            }
                                                        });
            }
            branchesRemoved = accumulate(threadBranchesRemoved.begin(), threadBranchesRemoved.end(), 0);
        }

        /*{
            // parallelly removing edges from graph. Need to lock nodes, because the same edges may be present to remove in many threads.
            WorkloadManager::parallelBlockExecution(0, Params::THREADS - 1, Params::THREADS, Params::THREADS,
                                                    [=, &edgesToRemove, &threadBranchesRemoved](unsigned a, unsigned b,
                                                                                                unsigned id) {
                                                        for (auto e : edgesToRemove[id]) {
                                                            G->lockNode(e.first);
                                                            if (G->removeDirectedEdge(e.first, e.second)) {
                                                                threadBranchesRemoved[id]++;
                                                            }
                                                            G->unlockNode(e.first);
                                                        }
                                                    });
            branchesRemoved = accumulate( threadBranchesRemoved.begin(), threadBranchesRemoved.end(), 0 );
        }*/

    }


    cerr << "\tdangling branches removed, " << branchesRemoved << " branches removed" << endl;

    GenomeStatisticsCollector::addData("dangling branches removed: ", branchesRemoved);

    for (int i = 0; i < G->size(); i++) (*G)[i].shrink_to_fit();

    TimeMeasurer::stopMeasurement("GraphSimplifier_removeDanglingBranches");
    return branchesRemoved;

}

int GraphSimplifier::removeDanglingBranchesFromNode(int beg, int maxOffset, vector<pair<int, int> > &edgesToRemove,
                                                    int thread_id) {

    int branchesRemoved = 0;

    vector<pair<int, int> > branchEnds;

    int infiniteLoopCheck = 0;
    VI neigh;

    unordered_map<int, int> dst, par;

    par[beg] = beg;
    for (int i = 0; i < (*G)[beg].size(); i++) {
        int v = (*G)[beg][i].first;
        par[v] = beg;

        was[thread_id][v] = true;
        neigh.push_back(v);

        int offset = G->getWeight(beg, i);

        infiniteLoopCheck = 0;

        while ((*G)[v].size() == 1) { // if v has only 1 son and was not visited yet

            int son = (*G)[v][0].first;
            if (was[thread_id][son]) break;

            was[thread_id][son] = true;
            neigh.push_back(son);

            par[son] = v;

            offset += G->getWeight(v, 0);

            if (offset < 0) {
                cerr << "offset < 0 in removing dangling braches" << endl;
                exit(1);
            }

            v = son;

            if (offset > maxOffset) break;

            if (infiniteLoopCheck++ > 1000000000) {
                cerr << "infinite loop in point 1" << endl;
                break;
            }
        }

        if ((*G)[v].size() == 0 && offset <= maxOffset) branchEnds.push_back({offset, v});

    }

    sort(branchEnds.begin(), branchEnds.end());

    int div = 0;
    if (branchEnds.size() == (*G)[beg].size()) div = 1;
    for (int i = 0; i < (int) branchEnds.size() - div; i++) {
        int v = branchEnds[i].second;

        infiniteLoopCheck = 0;
        while (v != beg) {

            branchesRemoved++;
            edgesToRemove.push_back({par[v], v});

            v = par[v];

            if (infiniteLoopCheck++ > 1000000000) {
                cerr << "infinite loop in point 2" << endl;
                break;
            }
        }
    }

    for (int a : neigh) {
        was[thread_id][a] = false;
    }
    par.clear();

    return branchesRemoved;
}


int GraphSimplifier::removeDanglingUpperBranches(int maxOffset) {
    G->reverseGraph();
    int branchesRemoved = removeDanglingBranches(maxOffset);
    G->reverseGraph();

    G->pruneGraph();

    TimeMeasurer::stopMeasurement("GraphSimplifier_removeDanglingUpperBranches");
    return branchesRemoved;
}


bool GraphSimplifier::contractPathNodes() {
    TimeMeasurer::startMeasurement("GraphSimplifier_contractPathNodes");

    bool anyContractionDone = false;
    int contractionsDone = 0;

    const bool USE_PARALLEL_CONTRACTION = true;

    if (!USE_PARALLEL_CONTRACTION) { // beginning of sequential path contraction

        VVPII GRev = G->getReverseGraphNeighborhoods();
        /**
         * Function does the same as removeDirected Edge in Graph class. It operates solely on the GRev VVPII structure. This is done to reduce memory peak.
         */
        auto removeDirectedEdge = [=, &GRev](int a, int b) {
            bool removed = false;
            int p = GRev[a].size() - 1;
            for (int i = GRev[a].size() - 1; i >= 0; i--) {
                if (GRev[a][i].first == b) {
                    swap(GRev[a][i], GRev[a][p]);

                    GRev[a].pop_back();
                    p--;
                    removed = true;
                }
            }

            return removed;
        };
        /**
        * Function does the same as removeDirected Edge in Graph class. It operates solely on the GRev VVPII structure. This is done to reduce memory peak.
        */
        auto addDirectedEdge = [=, &GRev](int a, int b, int offset) {
            if (a == b) return;
            VPII::iterator it = GRev[a].end();
            for (it = GRev[a].begin(); it != GRev[a].end(); ++it) if (it->first == b) break;
            if (it == GRev[a].end()) {
                GRev[a].push_back({b, offset});
            } else {
                int ind = it - GRev[a].begin();
                if (offset < GRev[a][ind].second) GRev[a][ind].second = offset;
            }
        };


        cerr << "Memory usage in contractPathNodes(), after creating reverse graph, before contracting nodes" << endl;
        MyUtils::process_mem_usage();

        int progressCounter = 0;
        deque<int> pathNodes;
        for (int i = 0; i < G->size(); i++) {
            if ((*G)[i].size() == 1 && GRev[i].size() == 1) pathNodes.push_back(i);
        }

        unsigned cnt = 0;
        while (!pathNodes.empty()) {
            int b = pathNodes.front();
            pathNodes.pop_front();

            if ((*G)[b].size() == 0) continue;

            int a = GRev[b][0].first;
            int c = (*G)[b][0].first;

            if (a == c) continue;


            if (G->contractPath(a, b, c)) {
                anyContractionDone = true;
                contractionsDone++;

                removeDirectedEdge(b, a);
                removeDirectedEdge(c, b);
                addDirectedEdge(c, a, G->findWeight(a, c));
            }

            if ((*G)[a].size() == 1 && GRev[a].size() == 1) {
                pathNodes.push_back(a);
            }
            if ((*G)[c].size() == 1 && GRev[c].size() == 1) {
                pathNodes.push_back(c);
            }

            cnt++;
            MyUtils::writeProgress(cnt, G->size(), progressCounter, "contracting paths progress",
                                   1); // just roughly accurate
        }
    } else { // parallel version
        const int THREADS_OLD = Params::THREADS;

        VB indeg1outdeg1(G->size(), true);
        VI *indegs = G->getInDegrees();
        for (unsigned i = 0; i < G->size(); i++) {
            if ((*indegs)[i] != 1) indeg1outdeg1[i] = false;
            if ((*G)[i].size() != 1) indeg1outdeg1[i] = false;
        }
        delete indegs;


        int blocks = 10 * Params::THREADS;
        VB threadAnyContractionDone(Params::THREADS, false);
        VI threadContractionsDone(Params::THREADS, 0);
        WorkloadManager::parallelBlockExecution(0, G->size() - 1, blocks, Params::THREADS,
                                                [=, &threadAnyContractionDone, &threadContractionsDone, &indeg1outdeg1](
                                                        unsigned x, unsigned y, unsigned id) {

                                                    for (unsigned i = x; i <= y; i++) {
                                                        if (indeg1outdeg1[i]) continue;

                                                        for (int j = 0; j < (*G)[i].size(); j++) {

                                                            int a = i;
                                                            PII p = (*G)[a][j];
                                                            int b = p.first;

                                                            if (!indeg1outdeg1[b]) continue;

                                                            int c = (*G)[b][0].first;
                                                            if (a == c) continue;

                                                            if (G->contractPath(a, b, c)) {
                                                                threadAnyContractionDone[id] = true;
                                                                threadContractionsDone[id]++;

                                                                j--; // need to check the same node as long as the contraction can be done
                                                            }
                                                        }
                                                    }

                                                });

        contractionsDone = accumulate(threadContractionsDone.begin(), threadContractionsDone.end(), 0);
        anyContractionDone = (contractionsDone > 0);

        Params::THREADS = THREADS_OLD; // #TEST
    }


    cerr << endl;
    cerr << "There were " << contractionsDone << " contractions" << endl;

    cerr << "Memory usage in contractPathNodes() after contracting nodes" << endl;
    MyUtils::process_mem_usage();

    TimeMeasurer::stopMeasurement("GraphSimplifier_contractPathNodes");
    return anyContractionDone;
}


int GraphSimplifier::mergeLength0Edges() {
    TimeMeasurer::startMeasurement("GraphSimplifier_mergeLength0Edges");
    Graph GRev = G->getReverseGraph();

    cerr << "Memory in mergeLength0Edges, after creating reverse graph" << endl;
    MyUtils::process_mem_usage();

    int nodesMerged = 0;
    int progressCounter = 0;
    for (int i = 0; i < G->size(); i++) {
        int a = i;

//        if( (*reads)[a] == nullptr ) continue;

        for (auto p : (*G)[i]) {
            int b = p.first;

//            if( (*reads)[b] == nullptr ) continue;

            if (p.second == 0 && G->containsEdgeShorterOrEqual(b, a, 0) &&
                /* THIS version here removes only nodes that have edge with offset 0 that forms a cycle of length 0 */
                //                ( Global::READS[a]->size() < Global::READS[b]->size() || ( Global::READS[a]->size() == Global::READS[b]->size() && a < b )  ) ){
                (Global::READS[a]->size() > Global::READS[b]->size() ||
                 (Global::READS[a]->size() == Global::READS[b]->size() && a < b))) {

                mergeNodes(a, b, GRev, 0);

                cerr << "merging nodes:" << endl << *(*reads)[a] << endl << *(*reads)[b] << endl
                     << "if done in GCPS, then it should be considered an error!" << endl;

                nodesMerged++;
                break;
            }
        }
        MyUtils::writeProgress(i + 1, G->size(), progressCounter, "mergeLEngth0Edges progress", 1);
    }

    cerr << "\tthere were " << nodesMerged << " nodes merged" << endl;
    TimeMeasurer::stopMeasurement("GraphSimplifier_mergeLength0Edges");

    return nodesMerged;
}

void GraphSimplifier::mergeNodes(int a, int b, Graph &GRev, int offsetAB) {

    for (PII x : (*G)[a]) GRev.removeDirectedEdge(x.first, a);

    G->mergeVertices(a, b, offsetAB);

    for (int j = 0; j < GRev[a].size(); j++) {
        int d = GRev[a][j].first;
        int w = GRev.getWeight(a, j);

        G->addDirectedEdge(d, b, w + offsetAB);

        GRev.addDirectedEdge(b, d, w + offsetAB);

        G->removeDirectedEdge(d, a);
    }

    G->clearNode(a);

    GRev.clearNode(a);
}

void GraphSimplifier::removeSmallOverlapEdges(int min_overlap_to_retain, int number_of_long_edges_to_retain) {
    vector<std::thread> parallelJobs;
    parallelJobs.reserve(Params::THREADS);

    int W = (int) ceil((double) G->size() / Params::THREADS);
    for (int i = 1; i < Params::THREADS; i++) {
        int a = min(i * W, G->size() - 1);
        int b = min((i + 1) * W - 1, G->size() - 1);

        parallelJobs.push_back(thread([=] {
            removeSmallOverlapEdgesJob(a, b, min_overlap_to_retain, number_of_long_edges_to_retain, i);
        }));
    }

    removeSmallOverlapEdgesJob(0, W - 1, min_overlap_to_retain, number_of_long_edges_to_retain, 0);

    for (auto &p : parallelJobs) p.join();
}

void GraphSimplifier::removeSmallOverlapEdgesJob(int a, int b, int min_overlap_to_retain,
                                                 int number_of_long_edges_to_retain, int thread_id) {
    VPII toRemove;
    set<int> zb;
    int progressCounter = 0;
    int removedEdges = 0;

    for (int i = a; i <= b; i++) {
        toRemove.clear();
        zb.clear();
        Read *r1 = Global::READS[i];

        for (auto &p : (*G)[i]) {
            Read *r2 = Global::READS[p.first];
            int overlap = Read::calculateReadOverlap(r1, r2, p.second);
            if (overlap < min_overlap_to_retain) toRemove.push_back({overlap, p.first});
        }

        sort(toRemove.begin(), toRemove.end());
        for (int r = 0; r < number_of_long_edges_to_retain && !toRemove.empty(); r++) toRemove.pop_back();

        for (auto &p : toRemove) zb.insert(p.second);
        for (int k = (*G)[i].size() - 1; k >= 0; k--) {
            int d = (*G)[i][k].first;
            if (zb.count(d)) {
                swap((*G)[i][k], (*G)[i].back());
                (*G)[i].pop_back();
                removedEdges++;
            }
        }

        (*G)[i].shrink_to_fit();

        if (thread_id == 0)
            MyUtils::writeProgress(i - a + 1, b - a + 1, progressCounter,
                                   "removing all except for at most " + to_string(number_of_long_edges_to_retain)
                                   + " edges with overlap  less than " + to_string(min_overlap_to_retain) +
                                   " in thread 0", 1);
    }

    if (thread_id == 0)
        cerr << endl << "There were " << removedEdges << " edges removed in thread " << thread_id << endl;
}

bool GraphSimplifier::removeEdgePairedEndSeparators() {

    cerr << "Proceeding to remove paired-end separators" << endl;
    VVPII GRev;
    {
        GRev = G->getReverseGraph().getV();
    }

    const int MAX_DST = 500;

    function<void(VVPII &, bool, int, int, unordered_set<int> &,
                  unordered_map<int, int> &was)> neighborhoodMarker = [=, &MAX_DST, &neighborhoodMarker]
            (VVPII &g, bool reversed, int v, int dst, unordered_set<int> &neigh, unordered_map<int, int> was) {


        if (dst > MAX_DST) return;


        if (was.count(v) && was[v] > dst) return;
        was[v] = dst;


        for (auto p : g[v]) {

            LPII edge;
            int tempDst = dst;
            if (reversed) {

                edge = G->getContractedEdgePath(p.first, v);

                for (PII x : edge) {
                    if (tempDst > MAX_DST) break;

                    neigh.insert(x.first);
                    tempDst += max(1, x.second);
                }
            } else {

                edge = G->getContractedEdgePath(v, p.first);
                VPII edge2(edge.begin(), edge.end());

                for (int j = (int) edge2.size() - 1; j >= 0 && tempDst < MAX_DST; j--) {
                    neigh.insert(edge2[j].first);
                    tempDst += max(1, edge2[j].second);
                }
            }

            if (tempDst < MAX_DST) neighborhoodMarker(g, reversed, p.first, tempDst, neigh, was);
        }
    };


    auto GV = G->getV();
    auto helper = [=, &GRev, &GV, &neighborhoodMarker](int a, int b) {

        unordered_set<int> neighB;
        unordered_map<int, int> was;
        neighborhoodMarker(GV, false, b, 0, neighB, was);

        was.clear();
        unordered_set<int> neighArev;
        neighborhoodMarker(GRev, true, a, 0, neighArev, was);

        int cnt = 0;
        const int THRESHOLD = 1;

        for (auto x : neighB) {
            if (neighArev.count(Read::getIdOfCompRevRead(Read::getIdOfPairedRead(x)))) {
                cnt++;
                if (cnt >= THRESHOLD) return false;
            }
        }

        return true;
    };


    const int SMALL_EDGE_LENGTH = 1.5 * Global::calculateAvgReadLength();
    VPII toRemove;

    auto helper2 = [=, &GRev, &SMALL_EDGE_LENGTH, &toRemove](int a, int b, int thread_id) {

        int progressCounter = 0;
        for (int i = a; i <= b; i++) {
//        for( int i=thread_id; i<G->size(); i+=Params::THREADS ){
            if ((*reads)[i] != nullptr) {

                for (PII d : (*G)[i]) {
                    if (d.second <= SMALL_EDGE_LENGTH) {
//                    if( d.second <= SMALL_EDGE_LENGTH && GRev[i].size() > 0 && (*G)[i].size() > 0 && ( GRev[i].size() >= 2 || (*G)[i].size() >= 2 ) ){

                        if (helper(i, d.first)) {

                            G->lockNode(0);
                            toRemove.push_back({i, d.first});
                            G->unlockNode(0);

                        }
                    }
                }
            }

            if (thread_id == 0)
                MyUtils::writeProgress(i - a + 1, b - a + 1, progressCounter, "removing paired-end separators", 1);
        }

    };


    vector<std::thread> parallelJobs;
    parallelJobs.reserve(Params::THREADS);

    int W = (int) ceil((double) G->size() / Params::THREADS);
    for (int i = 1; i < Params::THREADS; i++) {
        int a = min(i * W, G->size() - 1);
        int b = min((i + 1) * W - 1, G->size() - 1);

        parallelJobs.push_back(thread([=] { helper2(a, b, i); }));
    }

    helper2(0, W - 1, 0);

    for (auto &p : parallelJobs) p.join();

    cerr << endl << endl << "There are " << toRemove.size() << " edges to remove due to no paired-end neighborhood"
         << endl << endl;
    for (PII e : toRemove) {
        G->removeDirectedEdge(e.first, e.second);
    }


    if (!toRemove.empty()) return true;
    else return false;

}



