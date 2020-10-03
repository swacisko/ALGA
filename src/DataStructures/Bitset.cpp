/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Bitset.cpp
 * Author: sylwester
 * 
 * Created on November 20, 2018, 9:11 PM
 */

#include <DataStructures/Bitset.h>
#include <bitset>
#include <Utils/TimeMeasurer.h>
#include "Global.h"

#include "DataStructures/Bitset.h"

Bitset::Bitset(unsigned int size) : N(size) {
//    initializeStaticBlock();

//    if(size == 0) cerr << "size = 0 in bitset constructor" << endl;

//    blocks = BL_NUM(N-1) + 1;
//    V = VT( blocks() );
    if (blocks()) {
        V = new TYPE[blocks()];
        for (int i = 0; i < blocks(); i++) V[i] = ZEROS;
        lastElementModifier = ZEROS;
    }

//    for( long long i=0; i< ( IND_IN_BL(N)==0 ? BLOCK_SIZE : IND_IN_BL(N)  ); i++ ) lastElementModifier |= bits[i];
    for (long long i = 0; i < (IND_IN_BL(size) == 0 ? BLOCK_SIZE : IND_IN_BL(size)); i++)
        lastElementModifier |= bits.at(i);


}

Bitset::Bitset(vector<bool> &bits) : Bitset(bits.size()) {
    TYPE t = 0;
    unsigned ind = 0;
    for (unsigned i = 0; i < bits.size(); i++) {
        if (bits[i]) t |= ((TYPE) 1 << IND_IN_BL(i));

        if (IND_IN_BL(i + 1) == 0) {
            V[ind] = t;
            t = 0;
            ind++;
        }
    }

    if (t != ZEROS) {
        V[ind] = t;
    }
}

Bitset::Bitset(const Bitset &oth) {
//    blocks = oth.blocks;
    if (size() != oth.size()) {
        N = oth.N;
//        V.resize(oth.blocks());
        if (V != nullptr) {
            delete[] V;
            V = nullptr;
        }
        if (blocks()) V = new TYPE[oth.blocks()];
    }
    if (V == nullptr) V = new TYPE[blocks()];
    for (int i = 0; i < blocks(); i++) V[i] = oth.V[i];
//    V = oth.V;
    lastElementModifier = oth.lastElementModifier;
//    lastElementModifier = ZEROS;
//    for( long long i=0; i< ( IND_IN_BL(N)==0 ? BLOCK_SIZE : IND_IN_BL(N)  ); i++ ) lastElementModifier |= bits[i];
}

//
Bitset::~Bitset() {
    clear();
}

bool Bitset::operator!=(const Bitset &oth) {
    return !(*this == oth);
}

Bitset &Bitset::operator&=(const Bitset &oth) {
    int m = min(blocks(), oth.blocks());
    for (int i = 0; i < m; i++) V[i] &= oth.V[i];
    return *this;
}

Bitset Bitset::operator&(const Bitset &oth) {
    Bitset b(*this);
    b &= oth;
    return b;
}

/*
bool Bitset::operator<(const Bitset& oth) {
    
}

bool Bitset::operator<=(const Bitset& oth) {
    
}*/


Bitset Bitset::operator<<(int offset) {
    Bitset b(*this);
    b <<= offset;
//    return std::move(b);
    return b;
}

Bitset &Bitset::operator<<=(int offset) {
    try {
        int p = 0;
        const int beg = IND_IN_BL(offset);
        const int end = beg - 1;
        int q1 = BL_NUM(offset);
        int q2 = BL_NUM(offset - 1 + BLOCK_SIZE);

        while (q1 < blocks()) {
            TYPE a = (V[q1] & (zerosOnes[beg]));
//            TYPE a = (V.at(q1) & (zerosOnes.at(beg)));
            TYPE b = ZEROS;
            if (q2 < blocks()) b = (V[q2] & (~zerosOnes[end + 1]));
//            if (q2 < blocks()) b = (V.at(q2) & (~zerosOnes.at(end + 1)));


//            if( beg < 0 || beg >= BLOCK_SIZE || BLOCK_SIZE - 1 - end < 0 || BLOCK_SIZE-1-end >= BLOCK_SIZE  ){
////                cerr << "shifts: beg = " << beg << " BLOCKS_SIZE = " << BLOCK_SIZE << "   end = " << end << "    BLOCK_SIZE-1-end = " << BLOCK_SIZE-1-end << endl;
////                exit(1);
//            }

            if (beg < BLOCK_SIZE) a >>= beg;
            else a = 0;

            if (BLOCK_SIZE - 1 - end < BLOCK_SIZE) b <<= (BLOCK_SIZE - 1 - end);
            else b = 0;

            V[p] = (a | b);
//            V.at(p) = (a | b);

            p++;
            q1++;
            q2++;
        }

        while (p < blocks()) {
            V[p] = ZEROS;
//            V.at(p) = ZEROS;
            p++;
        }
        return *this;
    }
    catch (int e) {
        cerr << "EXCEPTION IN <<= in Bitset, caught e = " << e << endl;
        exit(1);
    }

}


bool Bitset::operator==(const Bitset &oth) {
    if (size() != oth.size()) return false;
    for (int i = 0; i < blocks(); i++) {
        if (V[i] != oth.V[i]) return false;
    }
    return true;
}

Bitset Bitset::operator>>(int offset) {
    Bitset b(*this);
    b >>= offset;
    return b;
}

Bitset &Bitset::operator>>=(int offset) {
    assert(blocks() >= 1);
    int p = blocks() - 1;

    const int end = BLOCK_SIZE - 1 - IND_IN_BL(offset);
    const int beg = end + 1;

    //  cerr << "beg: " << beg << "     end: " << end << endl;

    int q2 = p - BL_NUM(offset);
    int q1 = p - BL_NUM(offset - 1 + BLOCK_SIZE);

    while (q2 >= 0) {
        TYPE b = (V[q2] & (~zerosOnes[end + 1]));
//        TYPE b = ( V.at(q2) & ( ~zerosOnes.at(end+1) ) );
        TYPE a = ZEROS;
//        if( q1 >= 0 ) a = ( V.at(q1) & ( zerosOnes.at(beg) ) );
        if (q1 >= 0) a = (V[q1] & (zerosOnes.at(beg)));

        //   cerr << "p:     " << p << endl;
        //   cerr << "q1:    " << q1 << "      q2: " << q2 << endl << "a: " << toBinary(a) << endl << "b: " << toBinary(b) << endl << endl;
        //   cerr << "V[q1]: " << toBinary(V[q1]) << "   V[q2]: " << toBinary(V[q2]) << endl;


        if (beg < BLOCK_SIZE) a >>= beg;
        else a = 0;

        if (BLOCK_SIZE - 1 - end < BLOCK_SIZE) b <<= (BLOCK_SIZE - 1 - end);
        else b = 0;

        //   cerr << "a: " << toBinary(a) << endl << "b: " << toBinary(b) << endl << endl;

        V[p] = (a | b);
//        V.at(p) = (a|b);
        //   cerr << "V[p]: " << toBinary(V[p]) << endl;

        q1--;
        q2--;
        p--;
    }

    while (p >= 0) {
        V[p] = ZEROS;
        p--;
    }
//    while( p >= 0 ) V.at(p--) = ZEROS;

    return *this;
}

int Bitset::operator[](int ind) {
    if ((V[BL_NUM(ind)] & bits[IND_IN_BL(ind)]) != 0) return 1;
//    if( (V.at( BL_NUM(ind) ) & bits.at( IND_IN_BL(ind) ) ) != 0 ) return 1;
    else return 0;
}

Bitset &Bitset::operator=(const Bitset &oth) {
//    blocks = oth.blocks;
    if (size() != oth.size()) {
        N = oth.N;
//        if( oth.blocks() < 0 || oth.blocks() > 100 ){ cerr << "invalid in operator= in Bitset oth.blocks() = " << oth.blocks() << endl; exit(1); }
//        V.resize( (int)oth.blocks() );
        if (V != nullptr) delete[] V;
        V = new TYPE[blocks()];
    }
    if (V == nullptr) V = new TYPE[blocks()];
    for (int i = 0; i < blocks(); i++) V[i] = oth.V[i];
//    V = oth.V;

//    N = oth.N;
    lastElementModifier = oth.lastElementModifier;
//    lastElementModifier = ZEROS;
//    for( long long i=0; i< ( IND_IN_BL(N)==0 ? BLOCK_SIZE : IND_IN_BL(N)  ); i++ ) lastElementModifier |= bits[i];

    return *this;
}

Bitset &Bitset::set(int pos, bool value) {
    assert(blocks() >= 1);
    if (pos >= size()) return *this;

    if (value) V[BL_NUM(pos)] |= bits[IND_IN_BL(pos)];
//    if( value ) V.at( BL_NUM(pos) ) |= bits.at( IND_IN_BL(pos) );
    else V[BL_NUM(pos)] &= (~bits[IND_IN_BL(pos)]);
//    else V.at( BL_NUM(pos) ) &= (~bits.at( IND_IN_BL(pos) ) );

    if (BL_NUM(pos) == blocks() - 1) V[BL_NUM(pos)] &= lastElementModifier;
//    if( BL_NUM(pos) == blocks()-1 ) V.at( BL_NUM(pos) ) &= lastElementModifier;
    return *this;
}

Bitset &Bitset::set(int beg, int end, bool value) {
    assert(blocks() >= 1);
//    cerr << "before this = " << endl << *this << endl;
    if (end >= size()) end = max(beg, size() - 1);
    if (beg > end) return *this;

    int p = BL_NUM(beg);
    int q = BL_NUM(end);

//    cerr << "p = " << p << "   q = " << q << "   beg = " << beg << "    end = " << end <<  endl;

    TYPE A = V[p];
//    TYPE A = V.at(p);
    TYPE B = V[q];
//    TYPE B = V.at(q);
//    if( value ) A |= zerosOnes[ IND_IN_BL(beg) ];
    if (value) A |= zerosOnes.at(IND_IN_BL(beg));
//    else A &= (~zerosOnes[IND_IN_BL(beg)]);
    else A &= (~zerosOnes.at(IND_IN_BL(beg)));

    V[p] = A;
//    V.at(p) = A;

    for (int i = p + 1; i < q; i++) {
        if (value) V[i] = ONES;
//        if( value ) V.at(i) = ONES;
        else V[i] = ZEROS;
//        else V.at(i) = ZEROS;
    }


//    cerr << "p = " << p << "   q = " << q << endl;
//    cerr << "V[p]: " << endl << toBinary(V[p]) << endl;

//    if( value ) B |= ( IND_IN_BL(end+1) == 0 ? ~ZEROS :  (  ~zerosOnes[ IND_IN_BL(end+1) ] ) );
    if (value) B |= (IND_IN_BL(end + 1) == 0 ? ~ZEROS : (~zerosOnes.at(IND_IN_BL(end + 1))));
    else {
//        cerr << "p = " << p << "   q = " << q << "   B = " << endl << toBinary(B) << endl;
//        B &= ( IND_IN_BL(end+1) == 0 ? ZEROS : zerosOnes[ IND_IN_BL(end+1) ] );
        B &= (IND_IN_BL(end + 1) == 0 ? ZEROS : zerosOnes.at(IND_IN_BL(end + 1)));
//        cerr << "zerosOnes[end+1]:" << endl << toBinary(zerosOnes[ IND_IN_BL(end+1) ]) << endl;
//        cerr << "p = " << p << "   q = " << q << "   B = " << endl << toBinary(B) << endl;
    }


    if (p == q) {
        if (value) V[p] &= B;
//        if( value ) V.at(p) &= B;
        else V[p] |= B;
//        else V.at(p) |= B;
    } else V[q] = B;
//    else V.at(q) = B;

//    cerr << "V[p]: " << endl << toBinary(V[p]) << endl;

    if (q == blocks() - 1) V[q] &= lastElementModifier;
//    if( q == blocks()-1 ) V.at(q) &= lastElementModifier;

//    cerr << "after this = " << endl << *this << endl;
    return *this;
}

Bitset &Bitset::flip(int pos) {
    if (pos >= size()) return *this;

    if ((*this)[pos]) set(pos, false);
    else set(pos, true);
    return *this;
}

Bitset &Bitset::flip(int beg, int end) {
    assert(blocks() >= 1);
    if (end >= size()) end = max(beg, size() - 1);

    int p = BL_NUM(beg);
    int q = BL_NUM(end);

    int indBeg = IND_IN_BL(beg);
    int indEnd = IND_IN_BL(end);

    TYPE rest = V[p] & (~zerosOnes[indBeg]);
    //  cerr << "r: " << toBinary(rest) << endl;
    TYPE a = (V[p] & zerosOnes[indBeg]);
    //   cerr << "a: " << toBinary(a) << endl << endl;

    a |= (~zerosOnes[indBeg]);

    // cerr << "a: " << toBinary(a) << endl;
    a = ~a;
    // cerr << "a: " << toBinary(a) << endl;
    a |= rest;
    //   cerr << "a: " << toBinary(a) << endl;
    if (p == q) a &= (~zerosOnes[indEnd + 1]);
    //  cerr << "a: " << toBinary(a) << endl;

    for (int i = p + 1; i < q; i++) V[i] = ~V[i];

    rest = (V[q] & zerosOnes[indEnd + 1]);
    TYPE b = (V[q] & (~zerosOnes[indEnd + 1]));
    b |= zerosOnes[indEnd + 1];
    b = ~b;
    b |= rest;
    if (p == q) b &= (zerosOnes[indBeg]);

    //   cerr << "b: " << toBinary(b) << endl;

    if (p == q) V[p] = (a | b);
    else {
        V[p] = a;
        V[q] = b;
    }

    if (q == blocks() - 1) V[q] &= lastElementModifier;

    return *this;
}


Bitset Bitset::operator^(const Bitset &oth) {
    Bitset b(*this);
    b ^= oth;
    return b;
}

Bitset &Bitset::operator^=(const Bitset &oth) {
    int m = min(blocks(), oth.blocks());
    for (int i = 0; i < m; i++) V[i] ^= oth.V[i];
    return *this;
}

Bitset Bitset::operator|(const Bitset &oth) {
    Bitset b(*this);
    b |= oth;
    return b;
}

Bitset &Bitset::operator|=(const Bitset &oth) {
    int m = min(blocks(), oth.blocks());
    for (int i = 0; i < m; i++) V[i] |= oth.V[i];
    return *this;
}

Bitset Bitset::operator~() {
    assert(blocks() >= 1);
    Bitset b(*this);
//    for( auto &a : b.V ) a = ~a;
//    b.V.back() &= b.lastElementModifier;

    for (int i = 0; i < b.blocks(); i++) b.V[i] = ~b.V[i];
    if (b.size() > 0) b.V[b.blocks() - 1] &= lastElementModifier;

    return b;
}

int Bitset::count() const {
    int res = 0;
    if (sizeof(TYPE) == 4) {
//        for( auto a : V ) res += __builtin_popcount(a);
        for (int i = 0; i < blocks(); i++) res += __builtin_popcount(V[i]);
    } else if (sizeof(TYPE) == 8) {
//        for( auto a : V ) res += __builtin_popcountll(a);
        for (int i = 0; i < blocks(); i++) res += __builtin_popcountll(V[i]);
    }
    return res;
}


int Bitset::count(int a, int b) const {
    try {


        int res = 0;

        int p = BL_NUM(a);
        int q = BL_NUM(b);

        int beg = IND_IN_BL(a);
        int end = IND_IN_BL(b);


        if (BL_NUM(a) == BL_NUM(b)) {
            TYPE t = V[BL_NUM(a)];
//            TYPE t = V.at(BL_NUM(a));
            t &= zerosOnes[beg];
//            t &= zerosOnes.at(beg);
            t &= (IND_IN_BL(end + 1) == 0 ? ONES : ~zerosOnes[IND_IN_BL(end + 1)]);
//            t &= (IND_IN_BL(end + 1) == 0 ? ONES : ~zerosOnes.at(IND_IN_BL(end + 1)));
            if (sizeof(TYPE) == 4)res += __builtin_popcount(t);
            else if (sizeof(TYPE) == 8) res += __builtin_popcountll(t);
        } else {
            if (sizeof(TYPE) == 4)res += __builtin_popcount(V[p] & zerosOnes[beg]);
            else if (sizeof(TYPE) == 8) res += __builtin_popcountll(V[p] & zerosOnes[beg]);
//            else if (sizeof(TYPE) == 8) res += __builtin_popcountll(V.at(p) & zerosOnes.at(beg));

            for (int i = p + 1; i <= q - 1; i++) {
                if (sizeof(TYPE) == 4)res += __builtin_popcount(V[i]);
                else if (sizeof(TYPE) == 8) res += __builtin_popcountll(V[i]);
//                else if (sizeof(TYPE) == 8) res += __builtin_popcountll(V.at(i));
            }

            if (sizeof(TYPE) == 4)
                res += __builtin_popcount(V[q] & (IND_IN_BL(end + 1) == 0 ? ONES : ~zerosOnes[IND_IN_BL(end + 1)]));
            else if (sizeof(TYPE) == 8)
                res += __builtin_popcountll(V[q] & (IND_IN_BL(end + 1) == 0 ? ONES : ~zerosOnes[IND_IN_BL(end + 1)]));
//            else if (sizeof(TYPE) == 8)res += __builtin_popcountll(V.at(q) & (IND_IN_BL(end + 1) == 0 ? ONES : ~zerosOnes.at(IND_IN_BL(end + 1))));
        }
        return res;
    }
    catch (int e) {
        cerr << "caugh exception in count(a,b) in Bitset, e = " << e << endl;
        exit(1);
    }

}


int Bitset::size() const {
//    if(N<=0) cerr << "N <= 0: " << N << " in bitset, this should be impossible!" << endl;
    return N;
}

bool Bitset::any() {
//    for( auto a : V ){
    for (int i = 0; i < blocks(); i++) {
        if (V[i] != ZEROS) return true;
    }
    return false;
}

bool Bitset::all() {
//    for( auto a : V ){
    for (int i = 0; i < blocks(); i++) {
        if (V[i] != ONES) return false;
    }
    return true;
}


int Bitset::hash() const {
    long long P = 20000000001;
    long long res = 0;
//    for( auto a : V ){
    for (int i = 0; i < blocks(); i++) {
        auto a = V[i];
        res *= P;
        res += (a % P);
        if (res > P) res %= P;
    }
    return res;
}

VI Bitset::getAllSetBits() const {
    VI res;
    for (int i = 0; i < blocks(); i++) {
        TYPE b = V[i];
        while (b != 0) {
            TYPE t = (b & -b);
            int r;
            if (sizeof(TYPE) == 4) r = __builtin_ctz(t);
            else r = __builtin_ctzll(t);
            res.push_back(i * BLOCK_SIZE + r);
            b ^= t;
        }
    }
    return res;
}

void Bitset::clear() {
//    V.clear();
//    VT().swap(V);
    if (V != nullptr) delete[] V;
    V = nullptr;
}

int Bitset::lower_bound(int pos) {
    if (pos >= size()) return INF;
    TYPE b = V[BL_NUM(pos)];
    int d = IND_IN_BL(pos);

    if (d > 0) b &= zerosOnes[d];
    if (b != 0) {
        if (sizeof(TYPE) == 4) return (BL_NUM(pos) * BLOCK_SIZE + __builtin_ctz(b));
        else return (BL_NUM(pos) * BLOCK_SIZE + __builtin_ctzll(b));
    }

    for (int i = BL_NUM(pos) + 1; i < blocks(); i++) {
        if ((V[i] & ONES) != 0) {
            if (sizeof(TYPE) == 4) return i * BLOCK_SIZE + __builtin_ctz(V[i]);
            else return i * BLOCK_SIZE + __builtin_ctzll(V[i]);
        }
    }

    return INF;
}

unsigned int Bitset::upper_bound(unsigned int pos) {
    return lower_bound(pos + 1);
}

string Bitset::toBinary(TYPE a) {
    string s = "";
    while (a > 0) {
        if (a & 1 != 0) s += "1";
        else s += "0";
        a >>= 1;
    }
    while (s.size() < BLOCK_SIZE) s += "0";
    // reverse( s.begin(), s.end() );
    return s;
}


string Bitset::toString() {
    string s = "";
    for (int i = 0; i < size(); i++) {
        //  if( i > 0 && (i%BLOCK_SIZE==0) ) s += " ";
        if ((*this)[i] != 0) s += "1";
        else s += "0";
    }
    return s;
}

bool Bitset::operator<(const Bitset &oth) {
    assert(blocks() >= 1);
    int p = blocks() - 1;
    int q = oth.blocks() - 1;

    while (p > q) {
        if (V[p] > 0) return false;
        p--;
    }
    while (q > p) {
        if (oth.V[q] > 0) return true;
    }

    while (p >= 0) {
        if (V[p] != oth.V[p]) return V[p] < oth.V[p];
        p--;
    }
    return false;
}

ostream &operator<<(ostream &str, Bitset &b) {
    str << b.toString();
    return str;
}

void Bitset::negate() {
    assert(blocks() >= 1);
//    for( auto &a : V ){
    for (int i = 0; i < blocks(); i++) {
        auto b = V[i];
        V[i] = ~b;
    }
//    V.back() &= lastElementModifier;
    if (blocks()) V[blocks() - 1] &= lastElementModifier;
}


#define WRITE(myb, cppb){\
cerr << endl << myb << endl;\
for(int j=0; j<N; j++) cerr << cppb[j];\
cerr << endl;\
}

#define CHECK(myb, cppb) {\
\
TimeMeasurer::startMeasurement("myb");\
int mycnt = myb.count();\
TimeMeasurer::stopMeasurement("myb");\
\
TimeMeasurer::startMeasurement("cppb");\
int cppcnt = cppb.count();\
TimeMeasurer::stopMeasurement("cppb");\
bool same = true; \
for( int j=0; j<N; j++ ) if( myb[j] != cppb[j] ) same = false; \
if( same == false ){\
cerr << " diff after " << (i+1) << " cycles of operations" << endl;\
exit(1);\
}}

void Bitset::test() {

    const int N = 93;


    Bitset myb(N);
    bitset<N> cppb;

    Bitset temp(N);
    bitset<N> tempcpp;

    int operations = 30'000;
    for (int i = 0; i < operations; i++) {
        int a, b;

        Bitset temp2 = myb;
        myb = temp2;


//        swap(temp,myb);

        a = rand() % N;
        b = rand() % N;
        if (a > b) swap(a, b);

        // SET ON POSITION
        TimeMeasurer::startMeasurement("myb");
        myb.set(a, b % 2);
        TimeMeasurer::stopMeasurement("myb");

        TimeMeasurer::startMeasurement("cppb");
        cppb.set(a, b % 2);
        TimeMeasurer::stopMeasurement("cppb");

        //   cerr << "\rSET DONE" << flush;
        CHECK(myb, cppb)


        Bitset temp3(myb);
        myb = temp3;



        // SHL
        a = rand() % N;
        b = rand() % N;
        if (a > b) swap(a, b);

        TimeMeasurer::startMeasurement("myb");
//        myb <<= (b/10);
        myb <<= b;
        TimeMeasurer::stopMeasurement("myb");

        TimeMeasurer::startMeasurement("cppb");
//        cppb >>= (b/10);
        cppb >>= b;
        TimeMeasurer::stopMeasurement("cppb");

        //   cerr << "\rSHL DONE" << flush;
        CHECK(myb, cppb)




        // SET ON INTERVAL
        a = rand() % N;
        b = rand() % N;
        if (a > b) swap(a, b);

        TimeMeasurer::startMeasurement("myb");
        myb.set(a, b, b % 2);
        TimeMeasurer::stopMeasurement("myb");

        TimeMeasurer::startMeasurement("cppb");
        for (int k = a; k <= b; k++) cppb.set(k, b % 2);
        TimeMeasurer::stopMeasurement("cppb");

        //   cerr << "\rSET ON INTERVAL DONE" << flush;
        CHECK(myb, cppb)





        // FLIP ON INERVAL
        a = rand() % N;
        b = rand() % N;
        if (a > b) swap(a, b);

        TimeMeasurer::startMeasurement("myb");
        myb.flip(a, b);
        TimeMeasurer::stopMeasurement("myb");

        TimeMeasurer::startMeasurement("cppb");
        for (int k = a; k <= b; k++) cppb.flip(k);
        TimeMeasurer::stopMeasurement("cppb");

        //   cerr << "\rFLIP ON INTERVAL DONE" << flush;
        CHECK(myb, cppb)




        // SHR
        a = rand() % N;
        b = rand() % N;
        if (a > b) swap(a, b);

        TimeMeasurer::startMeasurement("myb");
        myb >>= (b / 10);
        TimeMeasurer::stopMeasurement("myb");

        TimeMeasurer::startMeasurement("cppb");
        cppb <<= (b / 10);
        TimeMeasurer::stopMeasurement("cppb");

        //    cerr << "\rSHR DONE" << flush;
        CHECK(myb, cppb)


        for (int k = 0; k < N; k++) {
            if (rand() % 2 == 1) {
                //   TimeMeasurer::startMeasurement("myb");
                temp.set(k, 1);
                //    TimeMeasurer::stopMeasurement("myb");

                //    TimeMeasurer::startMeasurement("cppb");
                tempcpp.set(k, 1);
                //    TimeMeasurer::stopMeasurement("cppb");
            } else {
                temp.set(k, 0);
                //    TimeMeasurer::stopMeasurement("myb");

                //    TimeMeasurer::startMeasurement("cppb");
                tempcpp.set(k, 0);
            }
        }



        // XOR AND OR SHIFT
        a = rand() % N;
        b = rand() % N;
        if (a > b) swap(a, b);

        TimeMeasurer::startMeasurement("myb");
        myb &= temp;
        myb ^= (temp << 1);
        myb |= (temp ^ (temp << 1));
        TimeMeasurer::stopMeasurement("myb");

        TimeMeasurer::startMeasurement("cppb");
        cppb &= tempcpp;
        cppb ^= (tempcpp >> 1);
        cppb |= (tempcpp ^ (tempcpp >> 1));
        TimeMeasurer::stopMeasurement("cppb");

        //   cerr << "\rXOR AND OR SHIFT DONE" << flush;
        CHECK(myb, cppb)



        // NEGATE
        a = rand() % N;
        b = rand() % N;
        if (a > b) swap(a, b);

        TimeMeasurer::startMeasurement("myb");
        myb = ~myb;
//        myb.negate();
        TimeMeasurer::stopMeasurement("myb");

        TimeMeasurer::startMeasurement("cppb");
        cppb = ~cppb;
        TimeMeasurer::stopMeasurement("cppb");

        CHECK(myb, cppb)


        a = rand() % N;
        b = rand() % N;
        if (a > b) swap(a, b);
        TimeMeasurer::startMeasurement("myb");
        int mycnt = myb.count(a, b);
        TimeMeasurer::stopMeasurement("myb");

        TimeMeasurer::startMeasurement("cppb");
        int cppcnt = 0;
        for (int k = a; k <= b; k++) if (cppb[k]) cppcnt++;
        TimeMeasurer::stopMeasurement("cppb");

        if (mycnt != cppcnt) {
            WRITE(myb, cppb);
            cerr << "cnt on interval diff in " << (i + 1) << " cycles" << endl;
            cerr << "mycnt = " << mycnt << "    cppcnt = " << cppcnt << endl;
            exit(1);
        }

    }

    cerr << "correct!" << endl;
    TimeMeasurer::writeAllMeasurements();
    exit(0);

}

int Bitset::mismatch(Bitset oth) {
    int m = min(blocks(), oth.blocks());
    for (int i = 0; i < m; i++) {
        if (V[i] != oth.V[i]) {
            if (sizeof(TYPE) == 4) return i * BLOCK_SIZE + __builtin_ctz(V[i]);
            else return i * BLOCK_SIZE + __builtin_ctzll(V[i]);
        }
    }
    return min(size(), oth.size());
}


Bitset::VT Bitset::bits;
Bitset::VT Bitset::zerosOnes;