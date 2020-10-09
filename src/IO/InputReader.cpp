/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InputReader.cpp
 * Author: sylwester
 * 
 * Created on November 23, 2018, 7:59 PM
 */

#include <Utils/TimeMeasurer.h>
#include <IO/InputReader.h>
#include <StatisticsGenerators/StatisticsGeneratorBigData.h>
#include <future>

#include "IO/InputReader.h"
#include "Global.h"
#include "StatisticsGenerators/GenomeStatisticsCollector.h"


string getComplimentaryString(string s) {
    string res = s;
    for (int i = 0; i < s.size(); i++) {
        if (Params::getNukl(s[i]) == Params::A) res[i] = 'T';
        if (Params::getNukl(s[i]) == Params::C) res[i] = 'G';
        if (Params::getNukl(s[i]) == Params::G) res[i] = 'C';
        if (Params::getNukl(s[i]) == Params::T) res[i] = 'A';
    }
    return res;

}


InputReader::InputReader() {
    STRreads = VI(Params::THREADS, 0);
    Nreads = VI(Params::THREADS, 0);

    NsInRead = VVI(Params::THREADS, VI(1000, 0));
}

//InputReader::InputReader(const InputReader& orig) {
//}

//InputReader::~InputReader() {
//}


void InputReader::readInput() {
    Params::inStream.open(Params::inStreamFilePath1);
    cin.rdbuf(Params::inStream.rdbuf());

    Global::READS.clear();
    cerr << "starting to read" << endl;
    readReads();
    cerr << "After first file there are " << Global::READS.size() << " reads" << endl;

    if (Params::ADD_PAIRED_READS && Params::INPUT_FILE_TYPE != Params::PFASTA) {
        if (readPairedReads()) {
            vector<Read *> reads(Global::READS.size());
            unsigned N = Global::READS.size();
            if (Params::ADD_COMP_REV_READS) {
                for (unsigned i = 0; 4ll * i < (LL) Global::READS.size(); i++) {
                    reads[i << 2] = Global::READS[i << 1];
                    reads[(i << 2) + 1] = Global::READS[(i << 1) + 1];
                    reads[(i << 2) + 2] = Global::READS[(N >> 1) + (i << 1)];
                    reads[(i << 2) + 3] = Global::READS[(N >> 1) + (i << 1) + 1];
                }
            } else {
                for (unsigned i = 0; 2ll * i < (LL) Global::READS.size(); i++) {
                    reads[i << 1] = Global::READS[i];
                    reads[(i << 1) + 1] = Global::READS[(N >> 1) + i];
                }
            }

            swap(Global::READS, reads);
            reads.clear();

            cerr << "After second file there are " << Global::READS.size() << " reads" << endl;
        }
    }

    for (unsigned i = 0; i < Global::READS.size(); i += 2) { // changing positions of read R and Rrevcomp
        swap(Global::READS[i], Global::READS[i + 1]); // will it swap anything?
    }


    for (unsigned i = 0; i < Global::READS.size(); i++) { // changing id to the new order.
        if (Global::READS[i] != nullptr) Global::READS[i]->setId(i);
    }

    Global::READS.shrink_to_fit();


    DEBUG(Global::READS.size());

    int l = Params::INF;
    int L = 0;
    double sumLength = 0;
    int cnt = 0;
    for (unsigned i = 0; i < Global::READS.size(); i++) { // some statistics
        if (Global::READS[i] == nullptr) continue;
        if (Global::READS[i]->size() == 0) {
            cerr << "read of size 0: " << endl << (*Global::READS[i]) << endl;
        }

        l = min(l, Global::READS[i]->size());
        L = max(L, Global::READS[i]->size());
        sumLength += Global::READS[i]->size();
        cnt++;
    }
    cerr << "sumLength: " << sumLength << endl;
    cerr << "max read length: " << L << endl;
    cerr << "min read length: " << l << endl;
    cerr << "avg read length: " << (sumLength / cnt) << endl;


    int Ns = 0, STRs = 0;
    for (int i = 0; i < Params::THREADS; i++) {
        Ns += Nreads[i];
        STRs += STRreads[i];
    }
    cerr << "There were " << 2 * Ns << " reads that contained N and were removed from graph creation process" << endl;
    cerr << "There were  " << 2 * STRs << " reads marked as STR and were removed from graph creation process" << endl;

    int rem = 0;
    for (int i = 0; i < Global::READS.size(); i++) if (Global::READS[i] == nullptr) rem++;
    cerr << "There are " << rem << " reads set to nullptr" << endl;

    cerr << endl << "Read with Ns:" << endl;
    for (int i = 0; i < 5; i++) {
        cnt = 0;
        for (int j = 0; j < Params::THREADS; j++) cnt += NsInRead[j][i];
        cerr << "There are " << 2 * cnt << " read with " << i << " N's" << endl;
    }

    cnt = 0;
    for (int j = 0; j < Params::THREADS; j++) cnt += NsInRead[j][349];
    cerr << "There are " << 2 * cnt << " read with " << 1 << " N's in the middle (at least 10 from end) of the read"
         << endl;

    GenomeStatisticsCollector::addData("different reads: ", Global::GRAPH.size());

}


string InputReader::readOneRead1(istream &str) {

    string s, empty;

    switch (Params::INPUT_FILE_TYPE) {
        case Params::MY_INPUT: {
            str >> s;
            break;
        }
        case Params::FASTA: {
            getline(str, empty);
            getline(str, s);
            break;
        }
        case Params::PFASTA: {
            getline(str, empty);
            getline(str, s);
            if (Params::ADD_PAIRED_READS == 0) {
                getline(str, empty);
                getline(str, empty);
            }

            break;
        }
        case Params::FASTQ: {
            getline(str, empty);
            getline(str, s);
            getline(str, empty);
            getline(str, empty);
            break;
        }
        default: {
            cerr << "DO NOT RECOGNIZE INPUT_FILE_TYPE" << endl;
        }

    }

    return s;
}


bool InputReader::readPairedReads() {
    Params::inStream.close();
    if (Params::inStreamFilePath2 == "") return false;
    cerr << "inStreamFilePath2 = " << Params::inStreamFilePath2 << endl;

    Params::inStream.open(Params::inStreamFilePath2);
    cin.rdbuf(Params::inStream.rdbuf());

    cerr << "reading paired reads" << endl;
    readReads();
    Params::inStream.close();

    return true;
}

void InputReader::readReads() {
    string s;


//    bool READ_MULTITHREAD = ( (Params::INPUT_FILE_TYPE == Params::MY_INPUT) ? false : true );
    bool READ_MULTITHREAD = true;


    if (READ_MULTITHREAD) {
        vector<vector<Read *> > parallelReads(Params::THREADS);
        vector<std::thread> parallelJobs;
        parallelJobs.reserve(Params::THREADS);
        for (int i = 1; i < Params::THREADS; i++) {
            parallelJobs.push_back(thread([=, &parallelReads] { readParallelJob(parallelReads, i); }));
        }

        readParallelJob(parallelReads, 0);
        for (auto &p : parallelJobs) p.join();

        cerr << "Parallel reading done, merging them to Global::READS" << endl;
//        for( int i=0; i<parallelReads.size(); i++ ) cerr << "There were " << parallelReads[i].size() << " reads read by " << i << "-th thread" << endl;

        int p = 0;
        while (true) {
            bool has = true;
            for (int i = 0; i < Params::THREADS; i++) {
                if (p < parallelReads[i].size()) {
                    if (p < 0 || p >= parallelReads[i].size()) {
                        cerr << "invalid p in merge reads to Global::READS" << endl;
                        exit(1);
                    }
                    Global::READS.push_back(parallelReads[i][p]);
                    if (Params::ADD_COMP_REV_READS) Global::READS.push_back(parallelReads[i][p + 1]);
                } else has = false;
            }
            if (Params::ADD_COMP_REV_READS) p += 2;
            else p++;

            if (!has) break;
        }


        cerr << "Reads merged to Global::READS" << endl;
    } else {
        int reads = 0;
        while ((s = readOneRead1(cin)) != "") { // this is for reading real data
            if (Params::REMOVE_READS_WITH_N && s.find("N") != string::npos) {
                int countN = 0;
                for (char c : s) if (c == 'N') countN++;
//                readsWithNSizes.push_back( countN );
                continue;
            }


            int p = 0;
            while (s[p] == ' ') p++;
            if (p > 0) s.erase(0, p);
            p = 0;
            while (p < s.size() && s[p] != ' ') p++;
            if (p < s.size()) s.erase(p, s.size() - p + 1);

            Global::addRead(s, Global::READS);
            if (Params::ADD_COMP_REV_READS)
                Global::addRead(MyUtils::getComplimentaryString(MyUtils::getReverse(s)), Global::READS);

            reads++;
            if (reads % 100000 == 0) cerr << "\rInputReader: " << reads << " already read" << flush;
        }
        cerr << endl;

    }

}

void InputReader::readParallelJob(vector<vector<Read *> > &reads, int thread_id) {

    string s;
    int readCount = 0;

    ifstream str;
    if (Global::READS.empty()) str.open(Params::inStreamFilePath1);
    else str.open(Params::inStreamFilePath2);

    for (int i = 0; i < thread_id; i++) readOneRead1(str);

    while ((s = readOneRead1(str)) != "") {

        int p = 0;
        while (s[p] == ' ') p++;
        if (p > 0) s.erase(0, p);
        p = 0;
        while (p < s.size() && s[p] != ' ') p++;
        if (p < s.size()) s.erase(p, s.size() - p + 1);

        if (thread_id < 0 || thread_id >= reads.size()) {
            cerr << "invalid thread_id in readParallelJobs" << endl;
            exit(1);
        }

        if (s.size() < Params::READ_END_TRIM_LEFT + Params::READ_END_TRIM_RIGHT + 10) {
//            cerr << "s = " << s << endl;
        } else {
            s.erase(s.begin(), s.begin() + min(Params::READ_END_TRIM_LEFT,
                                               (int) s.size())); // these two lines here were used, to produce best results!!! (porbably ends are not as good as they ought to be, even after trimmomatic).
            s.erase(s.end() - min(Params::READ_END_TRIM_RIGHT, (int) s.size()), s.end());
        }

        if (Params::TPN.find("rand_read_trim") != string::npos) {
            if (s[0] != 'A' && s[0] != 'C' && s[0] != 'G' && s[0] != 'T' && s[0] != 'N') {
                cerr << s << endl;
                exit(1);
            }
            int l = Params::getNuklNumber(Params::getNukl(s[0]));
            int r = Params::getNuklNumber(Params::getNukl(s.back()));
            s.erase(s.begin(), s.begin() + l);
            s.erase(s.end() - r, s.end());
        }


        bool containsN = false;

        int ncnt = 0;
        int nind = -1;
        for (int i = 0; i < s.size(); i++) {
            if (s[i] == 'N') ncnt++;

            if (s[i] != 'A' && s[i] != 'C' && s[i] != 'G' && s[i] != 'T' && s[i] != 'N' && s[i] != 'U') {
                cerr << "s[i] = " << s[i] << "   but should be A,C,G,T,N or U" << endl;
                exit(1);
            }
            if (s[i] == 'N' && Params::REMOVE_READS_WITH_N) {
                containsN = true;
                nind = i;
            }/* else if (s[i] == 'N') { // replace N'as with A'a
                s[i] = 'A';
            }*/
            else if (s[i] == 'N') {
                s[i] = Params::getNuklAsString(rand() & 3)[0];
//                ncnt++;
            } else if (Params::RNA && s[i] == 'U') s[i] = 'T';

        }

        NsInRead[thread_id][ncnt]++;
        if (ncnt == 1 && nind >= 10 && nind <= s.size() - 10) NsInRead[thread_id][349]++;

        int STR_THRESHOLD = 20;

        Read *r;
        if (Params::REMOVE_READS_WITH_N && containsN) {
            r = nullptr;
            Nreads[thread_id]++;
        } else {
            int threshold = Params::MIN_OVERLAP_AREA;
            threshold = STR_THRESHOLD;
            if (MyUtils::MinPeriod(s.c_str()) <= threshold) {
                r = nullptr;
                STRreads[thread_id]++;
            } else r = new Read(reads[thread_id].size(), s);
        }


        if (r != nullptr && r->size() == 0) {
            cerr << "read of size 0 corrersponding to s = " << s << endl;
        }

        reads[thread_id].push_back(r);

        if (Params::ADD_COMP_REV_READS) {
            reverse(s.begin(), s.end());
            s = getComplimentaryString(s);

            Read *r = nullptr;
            if (Params::REMOVE_READS_WITH_N && containsN) r = nullptr;
            else {
                int threshold = Params::MIN_OVERLAP_AREA;
                threshold = STR_THRESHOLD;
                if (MyUtils::MinPeriod(s.c_str()) <= threshold) r = nullptr;
                else r = new Read(reads[thread_id].size(), s);
            }

            reads[thread_id].push_back(r);
        }


        for (int i = 0; i < Params::THREADS - 1; i++) readOneRead1(str);

        readCount++;
        if (thread_id == 0) {
            if (readCount % 10'000 == 0) cerr << "\rInputReader: " << readCount << " already read by thread 0" << flush;
        }

    }

//    cerr << "thread " << thread_id << " read " << readCount << " reads" << endl;


    str.close();
}
