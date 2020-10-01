/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MyUtils.h
 * Author: sylwester
 *
 * Created on November 28, 2018, 1:59 PM
 */

#ifndef MYUTILS_H
#define MYUTILS_H

#include<iostream>
#include<sstream>
#include<vector>
#include <DataStructures/Graph.h>
#include <set>
#include <cstring>
#include "DataStructures/Read.h"

using namespace std;


#define WRITE1(A) { for(auto a : A) cerr << a << " "; cerr << endl;}
#define WRITE2(A) { for(auto a : A) { for(auto b : a) cerr << b << " "; cerr << endl; } }


class MyUtils {
public:
    MyUtils();

    MyUtils(const MyUtils &orig);

    virtual ~MyUtils();

    static int parseInt(string s);

    template<class _T>
    static string toString(_T a) {
        //  return to_string(a);
        stringstream ss;
        ss << a;
        return ss.str();
    }

    /**
     *
     * @param N
     * @return least X > N, such that X is prime
     */
    static int getNearestLowerPrime(int N);

    static string getComplimentaryString(string s);

    /**
     * Writes progress to cerr.
     * processed is the number of already processed elements
     * toProcess is the number of elements that are to be processed. If < 0 then i wil write processed every time that processed%accuracy == 0
     * progressCounter represents the number of percent already processed
     * text is the displayed text
     * accuracy: if 1 then i will write process state every 1 %. Otherwise i will write the state every 1/accuracy %.
     */
    static void
    writeProgress(unsigned processed, unsigned toProcess, int &progressCounter, string text, int accuracy = 1) {
        if (toProcess < 0) {
            if (processed % accuracy == 0) {
                cerr << "\r" << text << ": " << processed << flush;
            }
        } else {
            while (accuracy * 100 * (double) processed/*/(double)( toProcess )*/ >=
                   progressCounter * (double) (toProcess)) {
                if (accuracy == 1) cerr << "\r" << text << ": " << progressCounter << " %" << flush;
                else
                    cerr << "\r" << text << ": " << ((double) progressCounter / max((double) accuracy, 0.000001))
                         << " %" << flush;
                progressCounter++;
            }
        }
    }

    /*
     * replaces all occurences of given pattern in string str by string to
     */
    static string replaceAll(string str, const string &from, const string &to) {
        size_t start_pos = 0;
        while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
            str.replace(start_pos, from.length(), to);
            start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
        }
        return str;
    }

    static string getReverse(string s) {
        reverse(s.begin(), s.end());
        return s;
    }


    // Odwrotnosc modularna
// funkcja wyznacza takie x, ze a*x przystaje do 1 mod m
    static __int128 RevMod(__int128 a, __int128 m) {
        __int128 x, y;
        if (GCDW(a, m, x, y) != 1) return -1;
        // dokonaj przesuniecia zmiennej x, tak aby znajdowla sie w przedziale [0...m-1]
        x %= m;
        if (x < 0) x += m;
        return x;
    }


    /**
     * Creates and return complementary reverse contig to given one
     * @param ctg
     * @return
     */
    static Contig *getCompRevContig(Contig *ctg, vector<Read *> *reads) {

        vector<pair<Read *, int> > rcreads;

        auto V = ctg->getContainedReads();
        rcreads.push_back({(*reads)[V.back().first->getIdOfCompRevRead()], -1});

        for (int i = V.size() - 2; i >= 0; i--) {
            Read *r = V[i].first;
            Read *revcompR = (*reads)[r->getIdOfCompRevRead()];

            Read *r2 = V[i + 1].first;
            rcreads.push_back({revcompR, Read::getRightOffset(r, r2, V[i + 1].second)});
        }


        string s = MyUtils::getReverse(MyUtils::getComplimentaryString(ctg->getSequenceAsString()));
        Contig *c = new Contig(-1, s, rcreads);
//    c->correctSnipsInContig();

        return c;
    }


    static string getContractedPathAsString(vector<Read *> &reads, Graph &G, int a, int b) {
        LPII path = G.getContractedEdgePath(a, b);
        if (path.empty()) return "";
        else return getContractedPathAsString(reads, a, path);
    }

    //*******************************************************************************   MINIMALNY OKRES SLOWA
// funkcja jako parametr przyjmuje tekst, a zwraca dlugosc minimalnego okresu tego slowa
    static int MinPeriod(const char *str) { // DZIALA
        int pre[strlen(str) + 1], k = 0, q, m;
        pre[1] = 0;
        for (q = 1; str[q]; q++) {
            while (k > 0 && str[k] != str[q]) k = pre[k];
            if (str[k] == str[q]) k++;
            pre[q + 1] = k;
        }

        return strlen(str) - pre[strlen(str)];
    }

    // *************************************************************************** minimalna leksykograficzne cyklicznosc slowa
// funkcja jako rezultat zwraca numer pierwszej litery, rozpoczynajacej minimalna leksykograficznie cyklicznosc slowa
// ta funkcja to jest prawie to samo co znajdywanie maxymalnego sufiksu slowa, tylko zamiat "<"  mamy  ">"  i bierzemy modulo strlen(x)
    static int minLexCyc(const char *x) {
        int i = 0, j = 1, k = 1, p = 1, a, b, l = strlen(x);
        while (j + k <= (l << 1)) {
            if ((a = x[(i + k - 1) % l]) > (b = x[(j + k - 1) % l])) {
                i = j++;
                k = p = 1;
            } else if (a < b) {
                j += k;
                k = 1;
                p = j - i;
            } else if (a == b && k != p) k++;
            else {
                j += p;
                k = 1;
            }
        }
        return i;
    }


    // pre[k] oznacz dlugosc najdluzszego wlasciwego prefixo-sufixu slowa    t1-t2-t3-...-tk, gdzie t1-t2-...-tn to wzorzec
    // uwaga, pre jest liczone od 1 do wzo.size() a wzo to string, czyli od 0 do wzo.size()-1
    static vector<int> KMP(const char *str, const char *wzo) { // DZIALA
#define KMPH(z) while( k > 0 && wzo[k] != z[q] ) k = pre[k]; if( wzo[k] == z[q] ) k++;
        int pre[strlen(wzo) + 1], k = 0, q, m;
        pre[1] = 0;
        // wyznacz funkcje prefixowa dla zadanego wzorca
        for (q = 1; wzo[q]; q++) {
            KMPH(wzo);
            pre[q + 1] = k;
        }

        m = q;
        k = 0;

        VI indices;
        // dla kolejnych liter przeszukiwaneg o tekstu....
        for (q = 0; str[q]; q++) {
            // uzywajac funkcji prefixowej wyznacz dlugosc sufixu tekstu str[0...q], bedacego jednoczesnie sufixem wzorca
            KMPH(str);
            // jesli wyznaczony prefix jest dlugosci wzorca, to wywolujemy funkcje fun dla znalezionego wystapienia wzorca
            if (k == m) {
//                fun( q - m + 1 );
                indices.push_back(q - m + 1);
                k = pre[k];
            }
        }
        return indices;
    }

    static void process_mem_usage(double vm_usage = 0, double resident_set = 0);

private:

    static string getContractedPathAsString(vector<Read *> &reads, int beg, LPII &path) {
        string res = "";

        int head = beg;
        for (PII p : path) {
            for (int i = 0; i < p.second; i++) res += Params::getNuklAsString((*reads[head])[i]);
            head = p.first;
        }
        return res;
    }

// NWD euklidesem
    static __int128 GCD(__int128 x, __int128 y) {
        while (y > 0) swap(x %= y, y);
        return x;
    }

//funkcja wyznacza nwd oraz liczby l oraz k takie, ze l*a + k*b = nwd(a,b)
    static __int128 GCDW(__int128 a, __int128 b, __int128 &l, __int128 &k) {
        if (!a) {
            //gcd(0,b) = 0*0 + 1*b
            l = 0;
            k = 1;
            return b;
        }

        // Wyznacz rekurencyjnie wartosc najwiekszego wspolnego dzielnika oraz wspolczynniki
        // l oraz k

        __int128 d = GCDW(b % a, a, k, l);
//	cout << "Przed:  a = " << a << "   b = " << b << "   l = " << l << "   k = " << k << endl;
        //zaktualizuj wartosci wspolczynikow oraz zwróæ wynik
        l -= (b / a) * k;
//	cout << "Po:  b/a = " << b/a << "   a = " << a << "   b = " << b << "   l = " << l << "   k = " << k << endl;
        return d;
    }


};

#endif /* MYUTILS_H */

