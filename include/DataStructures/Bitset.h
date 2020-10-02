/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Bitset.h
 * Author: sylwester
 *
 * Created on November 20, 2018, 9:11 PM
 */

#ifndef BITSET_H
#define BITSET_H

#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
#include <cassert>

#define BL_NUM(index) ( (index) >> BLOCK_OFFSET )
// BL_NUM is the same as index / BLOCK_SIZE but faster

#define IND_IN_BL(index) ( (index) & MOD_MODIFIER )
// IND_IN_BL is just index % BLOCK_SIZE but faster

using namespace std;

typedef vector<int> VI;
typedef vector<long long> VLL;

class Bitset {
public:
//    typedef unsigned long long TYPE;
    typedef unsigned int TYPE;
    typedef vector<TYPE> VT;
    static const int BLOCK_SIZE = 8 * sizeof(TYPE);
    static const int BLOCK_OFFSET = (sizeof(TYPE) == 4) ? 5
                                                        : 6; // BLOCK_OFFSET is used to quickly calculate BL_NUM and IND_IN_BL
    static const int MOD_MODIFIER = ((1 << BLOCK_OFFSET) - 1);
    static VT bits; // bits[i] is the power 2^i. Bits[63] is -1 since this is 100000...
    static const TYPE ZEROS = (TYPE) 0;
    static const TYPE ONES = (~ZEROS);
    static const unsigned int INF = -1;
    static VT zerosOnes; // onesZeros[i] is the number that has i least significant bits set to 0, the rest are 1.


    Bitset(unsigned int size = 0); // creates a bitset with [size] bits accessible
    Bitset(const Bitset &orig);

    /**
     * Creates a bitset class for given set of bits.
     * @param bits
     */
    Bitset(vector<bool> &bits);

    ~Bitset();

    // void resize( long long s );


    Bitset operator^(const Bitset &oth);

    Bitset &operator^=(const Bitset &oth);

    Bitset operator&(const Bitset &oth);

    Bitset &operator&=(const Bitset &oth);

    Bitset operator|(const Bitset &oth);

    Bitset &operator|=(const Bitset &oth);

    Bitset operator<<(int offset); // translation to the left. after this V[0] = V[offset], and so on...
    Bitset &operator<<=(int offset);

    Bitset operator>>(int offset); // translation to the left. after this V[offset] = V[0], and so on...
    Bitset &operator>>=(int offset);

    bool operator<(
            const Bitset &oth); // CAUTION!! THIS OPERATOR TREATS BITSET AS BINARY NUMBER WITH REPRESENTATION INVERSE TO V. So V[0] is the lest significant bit;
    //   bool operator<=( const Bitset & oth );

    Bitset &operator=(const Bitset &oth);

    Bitset operator~();

    void negate();

    Bitset &flip(int pos);

    Bitset &flip(int beg, int end);

    Bitset &set(int pos, bool value); // sets given value at given position. Returns *this to enable chain operations
    Bitset &set(int beg, int end,
                bool value); // sets all bits in given range (both inclusive) to value. Returns *this to enable chain operations

    int operator[](int ind);

    bool operator==(const Bitset &oth);

    bool operator!=(const Bitset &oth);

    int hash() const; // returns value of this bitset (treated as number) modulo 10^9 + 1

    bool any();

    bool all();

    /**
     * Function returns the smallest index j such that *this and @{oth} differ at position j.
     * If one bitset is a 'prefix' of the other, then min(size(), oth.size()) is returned.
     * @param oth
     * @return
     */
    int mismatch(Bitset oth);

    int size() const;

    int count() const;

    int count(int a, int b) const; // counts set bits on given interval
    void clear();


    int lower_bound(
            int pos); // returns smallest index i >= pos such that i-th bit is 1. Works in O( d / BLOCK_SIZE ) where d is the distance between pos and first set bit after pos.
    unsigned int upper_bound(unsigned int pos);; // returns smallest index i > pos such that i-th bit is 1.


    VI getAllSetBits() const; // returns indexes of all bits that are 1. Works in O( (N/BLOCK_SIZE)  + #setBits )

    string toString();

    static string toBinary(TYPE a);

    friend ostream &operator<<(ostream &str, Bitset &b);


    static void test(); // tests whether random mass operations give the same results as c++ bitset


    static void initializeStaticBlock() {
        if (bits.empty()) {
            bits = VT(BLOCK_SIZE);
            bits[0] = 1;
            for (int i = 1; i < BLOCK_SIZE; i++) bits[i] = (bits[i - 1] << 1);

            zerosOnes = VT(BLOCK_SIZE + 1);
            zerosOnes[0] = ONES;
            for (int i = 1; i <= BLOCK_SIZE; i++) zerosOnes[i] = (zerosOnes[i - 1] << 1);
        }
    }

    TYPE getBlock(int b) { if (b >= blocks()) return ZEROS; else return V[b]; }

    long long countBlocks() { return blocks(); }

private:

//    VT V;
    TYPE *V = nullptr;
    TYPE lastElementModifier; // it is used to ensure that all bits in last element that have index > N are set to 0 after every operation
    unsigned N; // size of given vector
//    unsigned int N; // size of given vector
//    long long blocks;



    int blocks() const { return N == 0 ? 0 : BL_NUM(N - 1) + 1; }


};



/*
 
 __builtin_ffs(x)    this function returns 1 + lest significant 1-bit. If x == 0 then returns 0. e.g. for x=14 = (00..001110) returns 2;
 __builtin_ffsll(x)  for long long

 __builtin_clz(x)    returns number of leading 0-bits 
 __builtin_clzll(x)  for long long
 
 __builtin_ctz(x)    This function returns number of trailing 0-bits of x which starts from least significant bit position.
 __builtin_ctzll(x)  for long long

__builtin_popcount(x)   This function returns number of 1-bits of x.
__builtin_popcountll(x)   for long long
  
  
  
*/
#endif /* BITSET_H */

