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
#include <future>


ContigCreatorSinglePath::ContigCreatorSinglePath(Graph *G, vector<Read *> &reads) : ContigCreator(G, &reads) {
    GRev.clear();
}


vector<Contig *> ContigCreatorSinglePath::getAllContigs() {

    markReliablePredecessorsByPairedConnections();

    VB inDeg = G->hasPositiveIndegree();

    int progressCounter = 0;

    vector<Contig *> contigs;
    Contig::ID_COUNT = 0;


    { // parallel traversal
        auto getContigOmitShortCyclesFromJob = [=](unsigned a, unsigned b, unsigned thread_id) {
            vector<Contig *> contigs;
            for (unsigned i = a; i <= b; i++) {
                if ((*reads)[i] == nullptr) continue;
                if ((*G)[i].size() > 0) {
                    vector<Contig *> ctg = getContigOmitShortCyclesFrom(i);
                    contigs.insert(contigs.end(), ctg.begin(),
                                   ctg.end()); // if the vertex has indegree 0 and outdegree > 0 then i write it
                }
            }
            return contigs;
        };

        vector<std::future<vector<Contig *> > > futures(Params::THREADS - 1);

        VI nodesToCheck; // this is just the list of all nodes from which we will have to start creating paths. It is done to equally divide work into threads.
        for (unsigned i = 0; i < G->size(); i++) {
            if ((*reads)[i] == nullptr) continue;
            if ((*G)[i].size() > 0) nodesToCheck.push_back(i);
        }
        int WW = (int) ceil((double) nodesToCheck.size() / Params::THREADS);
        for (int i = 1; i < Params::THREADS; i++) {
            int a = nodesToCheck[i * WW];
            int b = min((i + 1) * WW - 1, (int) nodesToCheck.size() - 1);
            b = nodesToCheck[b];

            futures[i - 1] = std::async(std::launch::async, getContigOmitShortCyclesFromJob, a, b, i);
        }
        contigs = getContigOmitShortCyclesFromJob(nodesToCheck[0], nodesToCheck[WW - 1], 0);

        for (auto &p : futures) {
            auto thread_contigs = p.get();
            contigs.insert(contigs.end(), thread_contigs.begin(), thread_contigs.end());
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



    DEBUG(reliablePredecessors.size());
    unordered_map<int, unordered_set<int> >().swap(reliablePredecessors);
    unordered_map<unsigned, VPII>().swap(GRev);


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

    unordered_set<int> was;
    unordered_map<int, int> dst;
    dst[beg] = 0;

    was.insert(beg);
    VI visited(1, beg);

    vector<Contig *> contigs;
    vector<pair<Read *, int> > readsInContig;

    int predecessor = 0;
    for (int i = 0; i < (*G)[beg].size(); i++) {

        visited.clear();
        s = "";

        readsInContig.clear();
        readsInContig.emplace_back((*reads)[beg], -1);

        int offset = G->getWeight(beg, i);

        int p = (*G)[beg][i].first;
        predecessor = beg;

        addContractedPathToString(beg, p, s, readsInContig);

        was.insert(p);
        dst[p] = offset;
        visited.push_back(p);

        VPII nextCandidates = getNextStepCandidates(predecessor, p, readsInContig);
        int nextId = nextCandidates.empty() ? -1 : nextCandidates[0].first;
        offset = nextCandidates.empty() ? -1 : nextCandidates[0].second;
        int canBeNext = nextCandidates.size();

        if (canBeNext == 1) {
            if (offset != -1) dst[nextId] = dst[p] + offset;
            addContractedPathToString(p, nextId, s, readsInContig);
            predecessor = p;
            p = nextId;
        }

        while (canBeNext == 1) {
            visited.push_back(p);
            was.insert(p);

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

            if (p == -1 || was.count(p)) {
                break;
            }
        }

        assert(p >= 0 && p < reads->size());
        assert((*reads)[p] != nullptr);
        s += (*reads)[p]->getSequenceAsString();

        if (s.size() >= Params::CONTIG_MIN_OUTPUT_LENGTH) {
            G->lockNode(
                    0); // we aither need to lock some mutex or change type of static variable Contig::ID_COUNT to atomic type
            contigs.push_back(new Contig(Contig::ID_COUNT++, s, readsInContig));
            G->unlockNode(0);
        }

        if (canBeNext > 1) {
            if (!contigs.empty()) contigs.back()->setEndsInFork(true);
        }

        { // the same as clearing all nodes in visited
            was.clear();
            was.insert(beg);
            dst.clear();
            dst[beg] = 0;
        }
    }

    {
        visited.clear();
        s.clear();
    }

    return contigs;
}


VPII ContigCreatorSinglePath::getNextStepCandidates(int predecessor, int p, vector<pair<Read *, int>> &readsInContig) {

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

//    return false; // uncomment this line to FORBID EXTENDING PATHS
    if (reliablePredecessors.find(p) != reliablePredecessors.end() && reliablePredecessors[p].count(predecessor))
        return true;


//    if( was[d] == false || dst[p] - dst[ d ] + of > Params::CONTIG_CREATOR_SHORT_CYCLE_LENGTH ) return true;

    return false;
}


void
ContigCreatorSinglePath::addContractedPathToString(int a, int b, string &s, vector<pair<Read *, int> > &readsInContig) {
    LPII path = G->getContractedEdgePath(a, b);
    addContractedPathToString(a, path, s, readsInContig);
}

void
ContigCreatorSinglePath::addContractedPathToString(int a, LPII &path, string &s,
                                                   vector<pair<Read *, int>> &readsInContig) {

    for (auto p : path) {
        int offset = p.second;
        readsInContig.emplace_back((*reads)[p.first], p.second);

        assert((*reads)[a] != nullptr);
        for (int k = 0; k < offset; k++) {
            assert(k < (*reads)[a]->size());
            s += Params::getNuklAsString((*(*reads)[a])[k]);
        }
        a = p.first;
    }

}

void ContigCreatorSinglePath::markReliablePredecessorsByPairedConnections() {
    reliablePredecessors = unordered_map<int, unordered_set<int> >();


    minLengthOfEdgeForReliablePredecessor = Global::calculateAvgReadLength() * 2;


    {
        for (int i = 0; i < G->size(); i++) {
            for (PII p : (*G)[i]) {
                GRev[p.first].emplace_back(i, p.second);
            }
        }

        cerr << "Marking reliable predecessors, after creating GRev" << endl;
        MyUtils::process_mem_usage();
    }

    auto helper = [=](int a, int b, int thread_id) {
        for (int i = a; i <= b; i++) {

            auto it = GRev.find(i);
            if (it == GRev.end()) continue;

            if ((*G)[i].size() == 1 && (*G)[i][0].second >= minLengthOfEdgeForReliablePredecessor &&
                GRev[i].size() >= 1) {

                VI relPred = getReliablePredecessors(i);
                G->lockNode(0); // just locking so that many threads cannot access the same element at once
                if (!relPred.empty()) {
                    for (int x : relPred) reliablePredecessors[i].insert(x);
                }
                G->unlockNode(0);
            }
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
            cnt++;
        }

        s += to_string(DA[i].first) + " ";
    }

    s += "}";


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

