//
// Created by sylwester on 3/4/19.
//

#include <ContigCreators/ContigCreatorSinglePath.h>
#include <Utils/MyUtils.h>
#include <StatisticsGenerators/StatisticsGeneratorBigData.h>
#include <set>
#include <thread>
#include <Global.h>
#include <unordered_map>
#include <functional>

#include "ContigCreators/ContigCreatorSinglePath.h"

ContigCreatorSinglePath::ContigCreatorSinglePath(Graph *G, vector<Read *> &reads) : ContigCreator(G, &reads) {
    was = VB(G->size(), false);
    dst = VI(G->size(), 0);
}


vector<Contig *> ContigCreatorSinglePath::getAllContigs() {
    if (was.empty()) was = VB(G->size(), false);
    if (dst.empty()) dst = VI(G->size(), 0);

    markReliablePredecessorsByPairedConnections();

//    inDeg = G->getInDegrees();
    VB inDeg = G->hasPositiveIndegree();

    int progressCounter = 0;
//    int allSize = inDeg->size() + G->size();
    int allSize = 2 * G->size();

    vector<Contig *> contigs;
    Contig::ID_COUNT = 0;
//    for (int i = 0; i < inDeg->size(); i++) {
    for (int i = 0; i < G->size(); i++) {
        if ((*reads)[i] == nullptr) continue;

//        MyUtils::writeProgress(i + 1, inDeg->size() + G->size(), progressCounter, "creating contigs", 1);
        MyUtils::writeProgress(i + 1, G->size() + G->size(), progressCounter, "creating contigs", 1);

//        if ((*inDeg)[i] == 0 && (*G)[i].size() > 0) {
        if (inDeg[i] == false && (*G)[i].size() > 0) {
            vector<Contig *> ctg = getContigOmitShortCyclesFrom(i);

//            cerr << endl << "Found " << ctg.size() << " new contigs" << endl << endl;

            contigs.insert(contigs.end(), ctg.begin(),
                           ctg.end()); // if the vertex has indegree 0 and outdegree > 0 then i write it
//        } else if ((*inDeg)[i] == 0 && (*G)[i].size() == 0 && (*reads)[i]->size() >= Params::CONTIG_MIN_OUTPUT_LENGTH) {
        } else if (inDeg[i] == false && (*G)[i].size() == 0 &&
                   (*reads)[i]->size() >= Params::CONTIG_MIN_OUTPUT_LENGTH) {
            vector<pair<Read *, int> > containedReads;
            containedReads.push_back({(*reads)[i], 0});
            Contig *ctg = new Contig(Contig::ID_COUNT++, (*reads)[i]->getSequenceAsString(), containedReads);
            contigs.push_back(ctg); // if the vertex has indegree 0 and outdegree > 0 then i write it
        }
    }


    for (int i = 0; i < G->size(); i++) {
//        MyUtils::writeProgress(i + 1 + inDeg->size(), inDeg->size() + G->size(), progressCounter, "creating contigs", 1);
        MyUtils::writeProgress(i + 1 + G->size(), G->size() + G->size(), progressCounter, "creating contigs", 1);

        /*if( (*G)[i].size() >= 2  && (*inDeg)[i] > 0   ){
            vector<Contig*> ctg = getContigOmitShortCyclesFrom(i);
            contigs.insert( contigs.end(), ctg.begin(), ctg.end() );
        }
        else */
//        if ((*G)[i].size() >= 1 && (*inDeg)[i] > 0) {
        if ((*G)[i].size() >= 1 && inDeg[i]) {
            vector<Contig *> ctg = getContigOmitShortCyclesFrom(i);
            contigs.insert(contigs.end(), ctg.begin(), ctg.end());
        }
    }
    cerr << endl << "contigs created" << endl;
    cerr << "There are " << contigs.size() << " contigs after creation, before correcting SNPs" << endl;


    vector<std::thread> parallelJobs;
    parallelJobs.reserve(Params::THREADS);

    int W = (int) ceil((double) contigs.size() / Params::THREADS);
    for (int i = 1; i < Params::THREADS; i++) {
        int a = i * W;
        int b = min((i + 1) * W - 1, (int) contigs.size() - 1);
        parallelJobs.push_back(thread([=, &contigs] { correctSNPsJob(a, b, i, contigs); }));
    }
    correctSNPsJob(0, W - 1, 0, contigs);

    for (auto &p : parallelJobs) p.join();

    cerr << endl << "SNPs corrected" << endl;


//    delete inDeg;
//    inDeg = 0;


    DEBUG(reliablePredecessors.size());
    unordered_map<int, unordered_set<int> >().swap(reliablePredecessors);
//    VVPII().swap(GRev);
    unordered_map<unsigned, VPII>().swap(GRev);
    VI().swap(dst);
    VB().swap(was);


    cerr << "There are " << contigs.size() << " contigs returned" << endl;

    return contigs;
}

void ContigCreatorSinglePath::correctSNPsJob(int a, int b, int thread_id, vector<Contig *> &contigs) {

    int progressCounter = 0;
    for (int i = a; i <= b; i++) {
        Contig *ctg = contigs[i];

        ctg->correctSnipsInContig();
        if (thread_id == 0)
            MyUtils::writeProgress(i + 1 - a, b - a + 1, progressCounter, "correcting SNPs in contigs", 1);
    }
}

vector<Contig *> ContigCreatorSinglePath::getContigOmitShortCyclesFrom(int beg) {
    string s = "";

    int forks = 0;

//    unordered_map<int,bool> was;
//    unordered_map<int,int> dst;

    was[beg] = true;
    VI visited(1, beg);

    vector<Contig *> contigs;
    vector<pair<Read *, int> > readsInContig;

    int predecessor = 0;
    for (int i = 0; i < (*G)[beg].size(); i++) {

        visited.clear();
//        visited.push_back(beg);
        s = "";

        readsInContig.clear();
        readsInContig.emplace_back((*reads)[beg], -1);

        int offset = G->getWeight(beg, i);

        int p = (*G)[beg][i].first;
        predecessor = beg;

        addContractedPathToString(beg, p, s, readsInContig);


        was[p] = true;
        dst[p] = offset;
        visited.push_back(p);


        VPII nextCandidates = getNextStepCandidates(predecessor, p, readsInContig);
        int nextId = nextCandidates.empty() ? -1 : nextCandidates[0].first;
        offset = nextCandidates.empty() ? -1 : nextCandidates[0].second;
        int canBeNext = nextCandidates.size();

        if (canBeNext == 1) {
            if (offset != -1) dst[nextId] = dst[p] + offset;
            addContractedPathToString(p, nextId, s, readsInContig);
//            cerr << "appending path to contig " << Contig::ID_COUNT << endl;
            predecessor = p;
            p = nextId;
        }

        while (canBeNext == 1) {
            visited.push_back(p);
            was[p] = true;


            nextCandidates = getNextStepCandidates(predecessor, p, readsInContig);
            nextId = nextCandidates.empty() ? -1 : nextCandidates[0].first;
            offset = nextCandidates.empty() ? -1 : nextCandidates[0].second;
            canBeNext = nextCandidates.size();


            if (canBeNext == 1) {
                if (offset != -1) dst[nextId] = dst[p] + offset;
                addContractedPathToString(p, nextId, s, readsInContig);

                predecessor = p;
                p = nextId;
            }

            if (p == -1 || was[p]) {
                break;
            }

        }

        s += (*reads)[p]->getSequenceAsString();

        if (s.size() >= Params::CONTIG_MIN_OUTPUT_LENGTH) {
            contigs.push_back(new Contig(Contig::ID_COUNT++, s, readsInContig));
        }

        if (canBeNext > 1) {
            if (!contigs.empty()) contigs.back()->setEndsInFork(true);
        }


        for (int a : visited) {
            was[a] = false;
            dst[a] = 0;
        }

    }

    was[beg] = false;
    dst[beg] = 0;
    for (int a : visited) {
        was[a] = false;
        dst[a] = 0;
    }


    visited.clear();
    s.clear();

    return contigs;

}


VPII ContigCreatorSinglePath::getNextStepCandidates(int predecessor, int p, vector<pair<Read *, int>> &readsInContig) {
//    return getNextStepCandidatesByPairedReads( predecessor,p,readsInContig );

    VPII res;
    for (int j = 0; j < (*G)[p].size(); j++) {
        int w = (*G)[p][j].first;
        int of = (*G)[p][j].second;

        if (canBeNextStepCandidate(predecessor, p, w, of, readsInContig)) {
            res.push_back({w, of});
        }

    }

    return res;
}


bool ContigCreatorSinglePath::canBeNextStepCandidate(int predecessor, int p, int d, int of,
                                                     vector<pair<Read *, int>> &readsInContig) {

//        cerr << "\tCAUTION! do not lengthening paths!" << endl; return false; // #TEST

//    if( reliablePredecessors[p] == predecessor ) return true;
    if (reliablePredecessors[p].count(predecessor)) return true;





    /*bool REDUCE_MISASSEMBLIES = ( Params::TRAVERSE_TYPE == Params::TRAVERSE_SHALLOW );
    if( REDUCE_MISASSEMBLIES ){
        return false;
    }*/

//    if( was[d] == false || dst[p] - dst[ d ] + of > Params::CONTIG_CREATOR_SHORT_CYCLE_LENGTH ) return true;

    return false;

}


void
ContigCreatorSinglePath::addContractedPathToString(int a, int b, string &s, vector<pair<Read *, int> > &readsInContig) {
    LPII path = G->getContractedEdgePath(a, b);
//    if( path.size() > 100 ){
//        cerr << endl << "Found contracted edge path of size " << path.size() << endl;
//    }
    addContractedPathToString(a, path, s, readsInContig);
}

void
ContigCreatorSinglePath::addContractedPathToString(int a, LPII &path, string &s,
                                                   vector<pair<Read *, int>> &readsInContig) {

    for (auto p : path) {
        int offset = p.second;
        readsInContig.emplace_back((*reads)[p.first], p.second);

        for (int k = 0; k < offset; k++) s += Params::getNuklAsString((*(*reads)[a])[k]);
        a = p.first;
    }

}

void ContigCreatorSinglePath::markReliablePredecessorsByPairedConnections() {
//    reliablePredecessors = VI(G->size(),-1);
//    reliablePredecessors = vector<unordered_set<int>>(G->size());
    reliablePredecessors = unordered_map<int, unordered_set<int> >();


    minLengthOfEdgeForReliablePredecessor = Global::calculateAvgReadLength() * 2;
//    minLengthOfEdgeForReliablePredecessor = Global::calculateAvgReadLength();
//    minLengthOfEdgeForReliablePredecessor = 5;

    {
//        GRev = G->getReverseGraph().getV();
        for (int i = 0; i < G->size(); i++) {
            for (PII p : (*G)[i]) {
                GRev[p.first].emplace_back(i, p.second);
            }
        }

//        int cnt = 0;
//        for(int i=0; i<GRev.size(); i++) if( GRev[i].empty() ) cnt++;
//        DEBUG(cnt);
//        DEBUG(GRev.size());

        cerr << "Marking reliable predecessors, after creating GRev" << endl;
        MyUtils::process_mem_usage();
//        exit(1);
    }

    auto helper = [=](int a, int b, int thread_id) {
        int progressCounter = 0;
        for (int i = a; i <= b; i++) {

            auto it = GRev.find(i);
            if (it == GRev.end()) continue;

            if ((*G)[i].size() == 1 && (*G)[i][0].second >= minLengthOfEdgeForReliablePredecessor &&
                GRev[i].size() >= 1) {

//                int relPred = getReliablePredecessor(i);
//                G->lockNode(0); // just locking so that many threads cannot acces the same element at once
//                reliablePredecessors[i] = relPred;
//                if( relPred != -1 ){
//                    reliablePredecessors[i] = relPred;
//                    cerr << "Found reliable predecessor of node " << i << ": " << reliablePredecessors[i] << endl;
//                }
//                G->unlockNode(0);




                VI relPred = getReliablePredecessors(i);
//                VI relPred = VI( 1,getReliablePredecessor(i) ); // only one accepted predecessor, the same as above commented few lines
                G->lockNode(0); // just locking so that many threads cannot acces the same element at once
                if (!relPred.empty()) {
                    for (int x : relPred) reliablePredecessors[i].insert(x);
//                    cerr << "Found reliable predecessor of node " << i << ": " << reliablePredecessors[i] << endl;
                }
                G->unlockNode(0);

            }

//            MyUtils::writeProgress(i-a+1, b-a+1, progressCounter, "determination of reliable predecessors",1);
        }
    };


    vector<std::thread> parallelJobs;
    parallelJobs.reserve(Params::THREADS);

    int W = (int) ceil((double) G->size() / Params::THREADS);
    for (int i = 1; i < Params::THREADS; i++) {
        int a = i * W;
        int b = min((i + 1) * W - 1, (int) G->size() - 1);
        parallelJobs.push_back(thread([=] { helper(a, b, i); }));
    }

    helper(0, W - 1, 0);

    for (auto &p : parallelJobs) p.join();

}


int ContigCreatorSinglePath::getReliablePredecessor(int a) {

    int relPred = -1;
    int maxConn = 0;

    int b = (*G)[a][0].first; // successor of a

    for (PII pred : GRev[a]) {
        int d = pred.first; // predecessor of a
        int length = pred.second; // length of the edge (d,a)

        if (length < minLengthOfEdgeForReliablePredecessor) continue;

        int cnt = countPairedConnections(d, a, b);

        if (cnt > maxConn) {
            maxConn = cnt;
            relPred = d;
        }
    }


    if (maxConn < minPairedConnections) relPred = -1;

    return relPred;

}

VI ContigCreatorSinglePath::getReliablePredecessors(int a) {
    VI relPred;
    int maxConn = 0;

    int b = (*G)[a][0].first; // successor of a

    for (PII pred : GRev[a]) {
        int d = pred.first; // predecessor of a
        int length = pred.second; // length of the edge (d,a)

        if (length < minLengthOfEdgeForReliablePredecessor) continue;

        int cnt = countPairedConnections(d, a, b);

        if (cnt >= minPairedConnections) relPred.push_back(d);
    }

    return relPred;

}


int ContigCreatorSinglePath::countPairedConnections(int d, int a, int b) {

    auto edgeDA = G->getContractedEdgePath(d, a);
    auto edgeAB = G->getContractedEdgePath(a, b);

    // converting lists to vectors
    VPII DA(edgeDA.begin(), edgeDA.end());
    VPII AB(edgeAB.begin(), edgeAB.end());


    unordered_set<int> begOfAB;
    int dst = 0;

    string s = "\n\n\nbegOfAB: {";
    for (int i = 0; i < AB.size() && dst <= maxLengthOfInsertSize; i++) {
        dst += AB[i].second;
        begOfAB.insert(AB[i].first);

        s += to_string(AB[i].first) + " ";
    }

    s += "}\n\nEndOfDA: {";

    dst = 0;
    int cnt = 0;
    for (int i = (int) DA.size() - 1; i >= 0 && dst <= maxLengthOfInsertSize; i--) {
        dst += DA[i].second;

        int pairedId = Read::getIdOfPairedRead(DA[i].first);
        int pairedRevCompId = Read::getIdOfCompRevRead(pairedId);

        if (begOfAB.count(pairedId) || begOfAB.count(pairedRevCompId)) {
//            cerr << "Read:          " << *(*reads)[ DA[i].first ] << endl;
//            cerr << "Paired read:   " << *(*reads)[ pairedRevCompId ] << endl << endl;
            cnt++;
        }

        s += to_string(DA[i].first) + " ";
    }

    s += "}";

//    cerr << s << endl;

//    cerr << "cnt = " << cnt << endl;

    return cnt;

}


void ContigCreatorSinglePath::test() {

}

VPII ContigCreatorSinglePath::getNextStepCandidatesByPairedReads(int predecessor, int p,
                                                                 vector<pair<Read *, int>> &readsInContig) {
    VPII res;

    unordered_set<int> inContig;
    for (auto x : readsInContig) inContig.insert(x.first->getId());

    unordered_set<int> nodesForward;
    int dst = 0;
    int cnt = 0;


    function<bool(int, int, int, int, unordered_set<int> &)> helper = [=, &inContig, &helper](int p, int d, int dst,
                                                                                              int cnt,
                                                                                              unordered_set<int> &nodesForward) {

        cerr << "Entering helper p =  " << p << "   d = " << d << "   dst = " << dst << "   cnt = " << cnt << endl;
        if (cnt >= minPairedConnections) return true;
        if (dst > maxLengthOfInsertSize) return false;

        LPII edgePD = G->getContractedEdgePath(p, d);

        nodesForward.insert(d);


        for (auto x : edgePD) {
            int id = x.first;
            int pairedId = Read::getIdOfPairedRead(id);
            int pairedRevComp = Read::getIdOfCompRevRead(pairedId);
            int offset = x.second;

            if (inContig.count(pairedId) || inContig.count(pairedRevComp)) cnt++;
            dst += offset;

            if (cnt >= minPairedConnections) {
                nodesForward.erase(d);
                return true;
            }

            if (dst > maxLengthOfInsertSize) {
                nodesForward.erase(d);
                return false;
            }
        }

        for (PII neigh : (*G)[d]) {
            if (inContig.count(neigh.first) || nodesForward.count(neigh.first)) {
                nodesForward.erase(d);
                return false; // i do not want any cycles.
            }
        }

        for (PII neigh : (*G)[d]) {
            if (helper(d, neigh.first, dst, cnt, nodesForward)) return true;
        }

        return false;

    };


    for (PII neigh : (*G)[p]) {
        if (helper(p, neigh.first, dst, cnt, nodesForward)) {
            res.push_back(neigh);
        }
    }

    return res;

}

