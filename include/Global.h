/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Global.h
 * Author: sylwester
 *
 * Created on November 24, 2018, 6:03 PM
 */

#ifndef GLOBAL_H
#define GLOBAL_H

#include "DataStructures/Bitset.h"
#include "DataStructures/Graph.h"
#include "DataStructures/Read.h"
#include<cmath>
#include <unordered_set>
#include <unordered_map>


#define DEBUG(x)  cerr << #x << ": " << x << endl;

template<class _T, class _E>
ostream &operator<<(ostream &str, const pair<_T, _E> &pair) {
    str << "(" << pair.first << "," << pair.second << ")";
    return str;
}

template<class _T>
ostream &operator<<(ostream &str, const vector<_T> &vec) {
    int ile = 0;
    str << "{";
    for (auto p : vec) {
        if (ile > 0) str << ", ";
        ile++;
        str << p;
    }
    str << "}";
    return str;
}


template<class _T>
ostream &operator<<(ostream &str, const set<_T> &vec) {
    int ile = 0;
    str << "{";
    for (auto p : vec) {
        if (ile > 0) str << ", ";
        ile++;
        str << p;
    }
    str << "}";
    return str;
}


template<class _T>
ostream &operator<<(ostream &str, const unordered_set<_T> &vec) {
    int ile = 0;
    str << "{";
    for (auto p : vec) {
        if (ile > 0) str << ", ";
        ile++;
        str << p;
    }
    str << "}";
    return str;
}

template<class _T, class _E>
ostream &operator<<(ostream &str, const map<_T, _E> &vec) {
    int ile = 0;
    str << "{";
    for (auto p : vec) {
        if (ile > 0) str << ", ";
        ile++;
        pair<_T, _E> P = {p.first, p.second};
        str << P;
    }
    str << "}";
    return str;
}

template<class _T, class _E>
ostream &operator<<(ostream &str, const unordered_map<_T, _E> &vec) {
    int ile = 0;
    str << "{";
    for (auto p : vec) {
        if (ile > 0) str << ", ";
        ile++;
        pair<_T, _E> P = {p.first, p.second};
        str << P;
    }
    str << "}";
    return str;
}


struct pairhash {
public:
    template<typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const {
//        return std::hash<T>()(x.first) ^ (  std::hash<U>()(x.second) ^ (std::hash<T>( x.first )+1) ) ;
        return x.first ^ (x.second + 171);
    }
};


class Global {
public:


    Global();

    Global(const Global &orig);

    virtual ~Global();


    static int SEQUENCE_TOTAL_LENGTH;
    static vector<Read *> READS; // vector containing all reads.
    static Graph GRAPH; // here is the structure of created graph.

    static vector<char> pairedReadOffset;

    static void addRead(string s, vector<Read *> &READS = Global::READS);

    static double AVG_READ_LENGTH;

    static double calculateAvgReadLength() {
        AVG_READ_LENGTH = 0;
        int cnt = 0;
        for (Read *r : READS) {
            if (r != nullptr) {
                AVG_READ_LENGTH += r->size();
                cnt++;
            }
        }
        AVG_READ_LENGTH /= cnt;

        return AVG_READ_LENGTH;
    }


    /**
     * Function removes - that is deletes - read in Global::READS[id]
     * @param ids
     */
    static void removeRead(int id);

    /**
     * Function removes - deletes - all reads, that are represented by isolated vertices in graph GRAPH
     */
    static void removeIsolatedReads();

    static void generateFasta(string filename);

    static void writeReads(int a, int b); // writes all reads from given interval

    static void writeNodeWithNeighbors(int id) {
        cerr << id << endl;
        cerr << *READS[id] << endl;
        vector<pair<int, int>> neigh = GRAPH.getNeighbors(id);

        sort(neigh.begin(), neigh.end(), [](auto &&a, auto &&b) {
            if (a.second != b.second) return a.second < b.second;
            else return a.first < b.first;
        });
//        cerr << neigh[0].first << endl;

        for (int i = 0; i < neigh.size(); i++) {
//            cerr << neigh[i].first << "->" << neigh[i].second << endl;
            for (int k = 0; k < neigh[i].second; k++) {
                cerr << " ";
            }
            cerr << *READS[neigh[i].first] << endl;
        }
    }

    /*
     * Writes node dstTree with depth = 1
     */
    static void writeNodeTreeWithNeighbors(int id) {
        cerr << id << endl;
        cerr << *READS[id] << endl;
        vector<pair<int, int>> neigh = GRAPH.getNeighbors(id);

        sort(neigh.begin(), neigh.end(), [](auto &&a, auto &&b) {
            if (a.second != b.second) return a.second < b.second;
            elsereturn a.first < b.first;
        });

//        cerr << neigh[0].first << endl;

        for (int i = 0; i < neigh.size(); i++) {
            for (int k = 0; k < neigh[i].second; k++) {
                cerr << " ";
            }
            cerr << *READS[neigh[i].first] << endl;

            int b = neigh[i].first;

            vector<pair<int, int>> neigh2 = GRAPH.getNeighbors(b);
            sort(neigh2.begin(), neigh2.end(), [](auto &&a, auto &&b) {
                if (a.second != b.second) return a.second < b.second;
                elsereturn a.first < b.first;
            });

            for (int j = 0; j < neigh2.size(); j++) {
//                cerr << "    expanding" << endl;
                for (int k = 0; k < neigh2[j].second + neigh[i].second; k++) { cerr << " "; }
                cerr << *READS[neigh2[j].first] << endl;
            }


//            if( GRAPH[b].size() > 0 )cerr << "collapsing" << endl;

        }
    }


    static void writeNodeTreeWithNeighbors(int id, int depth, int offset = 0, int par = -1) {

        if (depth < 0) return;

        for (int i = 0; i < offset; i++) cerr << " ";
        cerr << *READS[id] << "    par = " << par << endl;

        vector<pair<int, int>> neigh = GRAPH.getNeighbors(id);
        sort(neigh.begin(), neigh.end(), [](auto &&a, auto &&b) {
            if (a.second != b.second) return a.second < b.second;
            elsereturn a.first < b.first;
        });

        for (int i = 0; i < neigh.size(); i++) {
            writeNodeTreeWithNeighbors(neigh[i].first, depth - 1, offset + neigh[i].second, id);
        }

    }


    static void writeReadPair(Read *r1, Read *r2, int offset) {
        cerr << *r1 << endl;
        for (int i = 0; i < offset; i++) cerr << " ";
        cerr << *r2 << endl;
    }

private:


};

#endif /* GLOBAL_H */

