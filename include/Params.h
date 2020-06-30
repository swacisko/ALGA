/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Params.h
 * Author: sylwester
 *
 * Created on November 23, 2018, 7:59 PM
 */

#ifndef PARAMS_H
#define PARAMS_H

#include "DataStructures/Bitset.h"
//#include "DataStructures/BigNum.h"
#include<fstream>
 // cannot include Global
#include<cmath>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

typedef vector<bool> VB;
typedef long long LL;
typedef vector<LL> VLL;

class Params {
public:
    
    Params();
    Params(const Params& orig);
    virtual ~Params();

    static const int INF = 1000000001;
    typedef long long KMER_HASH_TYPE;


    static int KMER_LENGTH_BUCKET;

    static int MOST_FREQUENTLY_USED_PARAMETER;
    static float SCALE;


    static int USE_LCS_LOW_ERROR_FILTER; // if 1 then i use AlgnmentControllerLowErrorRate and usual LCS. (it may be faster for low error tests). If 0 i use only usual, linear LCS.

    /*two reads will be marked as aligning if the LCS of their overlapping area has length at least MINIMAL_OVERLAP_RATE_FOR_LCS percent of overlapping area length*/
    static int MINIMAL_OVERLAP_RATE_FOR_LCS;

    // this is the threshold for AlignmentControllerLowErrorRate. If two read sequences have at least that percentage of overlap length in common, then they are considered as aligning.
    static int MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR; // ( (100 + MINIMAL_OVERLAP_RATE_FOR_LCS) >> 1 );

    /*
     * This is the maximal offset (int % of read length) that will be considered. If two reads will ahve offset greater than that value they will be deemed as not aligning
     * */
    static int MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT;

    /*
     * This is the maximal offset (in % of read length) for dangling branches. If a dangling branch has total offset (sum of lengths of edges) lower than this value, it will be removed from graph.
     */
    static int MAX_OFFSET_DANGLING_BRANCHES;


    /*
     * This is value in % of AVG_READ_LENGTH. So value ( MAX_OFFSET_PARALLEL_PATHS * AVG_READ_LENGTH / 100 ) will be maxOffset in removeParallelPaths in GraphSimplifier.
     */
    static int MAX_OFFSET_PARALLEL_PATHS;


    static int TRAVERSE_TYPE;
    static const int TRAVERSE_GREEDY = 0;
    static const int TRAVERSE_SHALLOW = 1;


    static int STR_REMOVAL_TYPE;
    static const int STR_MOA = 1;
    static const int STR_SIZE_HALF = 2;
    static const int STR_SIZE_MINUS_20 = 3;

    static int WRITE_STATISTICS; // if true (1) then i will collect and write statistics for given test. Otherwise only statistics that do not affect performance will be gathered and showed.

    /*
     * If true then every read in ninput is doubled - as the one in input and reverse complementary one.
     */
    static int ADD_COMP_REV_READS;
    static int ADD_REV_COMP_CONNECTIONS; // if true, then connections for revcomp edges will be added.
    static int ADD_PAIRED_READS;


    /*
        half of the size of the window during LCS check for each index. LCS will work PERFECTLY if there is an alignment in which two aligned
     *  nucleotides are at indices not more than MAX_ERROR_RATE_FOR_LCS apart.
     * E.g let us align two reads with offset 2:
     *  AAAACCCATTTTTG
     *    AAACGCTTTTGAA
     * resulted LCS is AAACCTTTTG
     * there is an alignment in which corresponding nucleotides in these two reads are not more than 1 index apart
     *  AAAC C TTTT G
     *   AAAC CTTTTG
     *
     */
    static int MAX_ERROR_RATE_FOR_LCS;

    /*
     * If true then all alignments will be tested only with AlignmentControllerLowErrorRate.
     * If false, then ACLER will be used only as positive filter.
     */
    static int USE_ACLER_INSTEAD_OF_ACLCS;

    /*
     * If 1 then i will use lexicographic index method (i will consider only kmers found by LI).
     */
    static const int USE_LI = 1;

    /*
     * If true, then i will create graph using pairwise-kmer method. Othwerwise i will use usual alignment controller (bitmap + LCS).
     */
    static int PAIRWISE_KMER_BRANCH_GRAPH_CREATION;
    static int PREF_SUF_GRAPH_CREATION;
    static int ALGORITHM_IN_USE;

    static int USE_GRAPH_CREATOR_SUPPLEMENT;

    static int REMOVE_SHORT_PARALLEL_PATHS_MST_FAST;

    /*
     * This is the minimal overlap (in nucleotides, not percent of length) that two reads must have to be counted as aligning.
     */
    static int MIN_OVERLAP_PREF_SUF;

    /*
     * This is the minimal offset that two reads must have to be aligned. If this is 0 then reaads with offset 0 will be aligned. Otherwise no cycles of length 0 will occur.
     */
    static int MIN_OFFSET_FOR_ALIGNMENT;

    /*
     * This is the threshold for accepting a pair of reads with at least MIN_PAIRWISE_KMER_OVERLAP % of total number of its kmers as overlapping. It will be used in GraphCreatorPairwiseKmer.
     */
  //  static const int MIN_PAIRWISE_KMER_OVERLAP = 5;

    /*
     * This is the upper bound for a value of a hash. If a hash is long and would take larger values, they will be taken modulo this number
     */
    static KMER_HASH_TYPE MAX_HASH_CONSIDERED;

    /**
     * MIN_OVERLAP is equal to the minimum overlap are of 2 reads to consider them as aligning.
     */
    static int MIN_OVERLAP_AREA;

    /**
     * TestParameterName - string representing name of the variable tested
     */
    static string TPN;

    /**
     * TestParameterValue - used for testing particualr parameters
     */
    static string TPV;


    static int LI_KMER_LENGTH;
    static int LI_KMER_INTERVALS;
//    static vector<int> nucleotidePriorities; // nucleotidePriorities[i] is the number that will have i-th type of nucleotide. This is used to create different


    static int GRAPH_CREATOR_PREF_SUF_2_REMOVAL_FREQUENCY;

    static int CONTIG_MIN_OUTPUT_LENGTH; // if created contig will be shorter than this value it will not be written and will not be used in statistics. (what for if reads are as long as these contigs? )
    static int SWAT_BRANCH_WINDOW_SIZE; // size of the window in alignment controller. this is the depth of checking into each branch.


    static int USE_STREAMS_INPUT;
    static int REMOVE_READS_WITH_N;
    static int ALLOW_INCLUDED_ALIGNMENT; // if true then shorter reads that are included within longer ones (as a k-mer; with possible errors of course) will be ignored

    static bool REDIRECT_CERR;

    static int REASSEMBLE_CONTIGS;

    static int VISUALIZE_CONTIGS;
    static int SCAFFOLD_CONTIGS;

    /**
     * This percentage of reads will be retained. This may be usefull if there are huge data sets with very large cover.
     */
    static int RANDOM_READS_TO_RETAIN_PERCENTAGE;

    /*
     * This is value that will be used in mergeLength2Cycles in GraphSimplifier.
     */
    static int CYCLES_MERGE_THRESHOLD;

    /*
     * If true, than truly metric triangles will also be removed.
     */
    static int REMOVE_TRULY_METRIC_TRIANGLES;

    /*
     * set this to INF if do not want to use maps in graph at all. This may be essential if we wanto to access graph state during graph creation.
     * In GraphCreatorSwatBranch and GraphCreatorPairwiseKmer we fortuantely do not need to access graph state durign construction - we only add edges.
    */
    static int GRAPH_NONMAP_SIZE_THRESHOLD;

    static int GRAPH_SIMPLIFIER_REMOVE_ITERATIONS;

    static int USE_SERIALIZATION_FOR_BUCKETS;
    static int SERIALIZE_GRAPH_BEFORE_SIMPLIFIER; // graph wull be serialized before simplifier
    static int SERIALIZE_GRAPH_AFTER_METRIC_TRIANGLES; // graph will be serialized after simplifier
    static int SERIALIZE_GRAPH_AFTER_SIMPLIFIER; // graph will be serialized after simplifier
    static int DESERIALIZE_GRAPH; // if true then alignment graph will be deserialized if available.
    static int WRITE_IN_GRASSHOPPER_FORMAT;
  //  static const int USE_ALREADY_SERIALIZED = 0;
    static int REMOVE_SERIALIZED_FILES_AFTER_USE; // removes serialized bucket with kmers after using them
    static string TEST_NAME;
    static string RUN_DATE;

    static const int MY_INPUT = 1;
    static const int FASTA = 2;
    static const int FASTQ = 3;
    static const int PFASTA = 4;


    static int RNA;

    static int REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP;
    static int REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN;

    static int GENERATE_FASTA;
    static int GENERATE_PFASTA;

    static int READ_END_TRIM_LEFT;
    static int READ_END_TRIM_RIGHT;

    /*
     * If 1, then reads will be corrected and written in seperate fasta files, and assembly performed. If 2, then no assembly will be done.
     */
    static int CORRECT_READS;


    static int THREADS;
    static int MAX_OFFSET_FOR_RETAINING_PAIRED_END_EDGE;
    static int LI_PRIORITIES_TO_CONSIDER; // max is 4. This number of priorities will be considered in order A,C,G,T.

    /**
     * If considered contig contains S reads, and NEW_READS_PER_CONTIG_PERCENTAGE/100 * S of them were not contained in any other contig written earlier,
     * then this contig will considered as proper contig. Otherwise it will be considered as a revcomp of a contig considered before.
     */
    static int NEW_READS_PER_CONTIG_PERCENTAGE;

    static int INPUT_FILE_TYPE;

    static int CONTIG_CREATOR_SHORT_CYCLE_LENGTH;

    static int JUST_CREATE_ALIGNMENT_GRAPH_AND_BASIC_STATISTICS; // if true, then i will not write any contigs nor do all simplifications in GraphSimplifier.
    static int COMPARE_GRAPHS;

    static int USE_DATE_IN_STREAM_FILE_NAMES; // if true then stdout and stderr will be redirected to a file with name that contains run date.

    static int MIN_OVERLAP_RATE;

    static int ALIGNMENT_CONTROLLER_SAME_ENDS_LENGTH;

    static const int PREF_READS_ONLY_DUPLICATES = 1;
    static const int PREF_READS_ALL_PREFIX_READS = 2;
    static const int PREF_READS_NONE = 3;
    static int REMOVE_PREF_READS_TYPE;

    typedef int NUKL_TYPE;
    static const NUKL_TYPE A = 0;
    static const NUKL_TYPE C = 1;
    static const NUKL_TYPE G = 2;
    static const NUKL_TYPE T = 3;
    static const int N = 15;

    static string getNuklAsString( int nukl );

    static int getNukl( char a ); // for given character returns int corresponding to the nucleotide
    static int getNuklNumber( NUKL_TYPE nucl); // see definition - it will be clear.

    static void writeParams();
    static void initializeParams( int argc, char** argv );

    static void createOutputFileName();
    

    static void closeStreams(){ inStream.close(); outStream.close(); errStream.close(); }

    static ifstream inStream;
    static ofstream outStream;
    static ofstream errStream;
    static string fileExtension;
    static string inStreamFilePath1;
    static string inStreamFilePath2;

    static string outStreamFileName;

private:



};

#endif /* PARAMS_H */

