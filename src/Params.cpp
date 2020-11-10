/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Params.cpp
 * Author: sylwester
 * 
 * Created on November 23, 2018, 7:59 PM
 */

#include <Params.h>
#include <Utils/MyUtils.h>
#include <chrono>

#include "Params.h"
#include "DataStructures/Kmer.h"

Params::Params() {
}

Params::Params(const Params &orig) {
}

Params::~Params() {
}

void Params::writeParams() {
    cerr << endl << "PARAMS:" << endl;

    cerr << "LI_KMER_LENGTH: " << LI_KMER_LENGTH << endl;
    cerr << "LI_KMER_INTERVALS: " << LI_KMER_INTERVALS << endl;
    cerr << "MIN_OVERLAP_PREF_SUF: " << MIN_OVERLAP_PREF_SUF << endl;

    cerr << "KMER_LENGTH_BUCKETS: " << KMER_LENGTH_BUCKET << endl;

    cerr << "MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT: " << MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT << endl;
    cerr << "MAX_ERROR_RATE_FOR_LCS: " << MAX_ERROR_RATE_FOR_LCS << endl;

    cerr << "USE_LCS_LOW_ERROR_FILTER: " << USE_LCS_LOW_ERROR_FILTER << endl;
    cerr << "MINIMAL_OVERLAP_RATE_FOR_LCS: " << MINIMAL_OVERLAP_RATE_FOR_LCS << endl;
    cerr << "MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR: " << MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR << endl;
    cerr << "USE_ACLER_INSTEAD_OF_ACLCS: " << USE_ACLER_INSTEAD_OF_ACLCS << endl;

    cerr << "USE_LI: " << USE_LI << endl;


    cerr << "MAX_OFFSET_PARALLEL_PATHS: " << MAX_OFFSET_PARALLEL_PATHS << endl;
    cerr << "MAX_OFFSET_DANGLING_BRANCHES: " << MAX_OFFSET_DANGLING_BRANCHES << endl;


    cerr << "ADD_COMP_REV_READS: " << ADD_COMP_REV_READS << endl;
    cerr << "ADD_PAIRED_READS: " << ADD_PAIRED_READS << endl;
    cerr << "INPUT_TYPE: "
         << (INPUT_FILE_TYPE == MY_INPUT ? "MY_INPUT" : (INPUT_FILE_TYPE == FASTA ? "FASTA" : "FASTQ")) << endl;
    cerr << "(LL)MAX_HASH_CONSIDERED: " << (LL) MAX_HASH_CONSIDERED << endl;
    cerr << "CONTIG_MIN_OUTPUT_LENGTH: " << CONTIG_MIN_OUTPUT_LENGTH << endl;
    cerr << "CONTIG_CREATOR_SHORT_CYCLE_LENGTH: " << CONTIG_CREATOR_SHORT_CYCLE_LENGTH << endl;

    cerr << "WRITE_STATISTICS: " << WRITE_STATISTICS << endl;

    cerr << "SERIALIZE_GRAPH_BEFORE_SIMPLIFIER: " << SERIALIZE_GRAPH_BEFORE_SIMPLIFIER << endl;
    cerr << "SERIALIZE_GRAPH_AFTER_SIMPLIFIER: " << SERIALIZE_GRAPH_AFTER_SIMPLIFIER << endl;
    cerr << "DESERIALIZE_GRAPH: " << DESERIALIZE_GRAPH << endl;

    cerr << "NEW_READS_PER_CONTIG_PERCENTAGE: " << NEW_READS_PER_CONTIG_PERCENTAGE << endl;
    cerr << "ALIGNMENT_CONTROLLER_SAME_ENDS_LENGTH: " << ALIGNMENT_CONTROLLER_SAME_ENDS_LENGTH << endl;

    cerr << "REMOVE_READS_WITH_N: " << REMOVE_READS_WITH_N << endl;

    cerr << "MIN_OFFSET_FOR_ALIGNMENT: " << MIN_OFFSET_FOR_ALIGNMENT << endl;


    cerr << "REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN: " << REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN << endl;


    cerr << "A = " << A << "  C = " << C << "   G = " << G << "   T = " << T << "   N = " << N << endl;
    cerr << "RNA: " << RNA << endl;

    cerr << "TRAVERSE_TYPE: " << TRAVERSE_TYPE << endl;


    cerr << "TEST_NAME: " << TEST_NAME << endl;
    cerr << "RUN_DATE: " << RUN_DATE << endl;

    cerr << endl << endl;

    cerr << "MOST_FREQUENTLY_USED_PARAMETER: " << MOST_FREQUENTLY_USED_PARAMETER << endl;
    cerr << "REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP: " << REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP << endl;
    cerr << "READ_END_TRIM_LEFT: " << READ_END_TRIM_LEFT << endl;
    cerr << "READ_END_TRIM_RIGHT: " << READ_END_TRIM_RIGHT << endl;

    cerr << "PREF_READS_ONLY_DUPLICATES: " << PREF_READS_ONLY_DUPLICATES << endl;
    cerr << "PREF_READS_ALL_PREFIX_READS: " << PREF_READS_ALL_PREFIX_READS << endl;
    cerr << "PREF_READS_NONE: " << PREF_READS_NONE << endl;
    cerr << "REMOVE_PREF_READS_TYPE: " << REMOVE_PREF_READS_TYPE << endl;
    cerr << "ALGORITHM_IN_USE: "
         << (ALGORITHM_IN_USE == PREF_SUF_GRAPH_CREATION ? "PS-no errors admitted" : "PKB-erros admitted") << endl;
    cerr << "THREADS: " << THREADS << endl;

    cerr << endl;
    cerr << "Test Parameter Name: " << TPN << endl;
    cerr << "Test Parameter Value: " << TPV << endl;

    cerr << endl << endl;
}

int Params::getNuklNumber(Params::NUKL_TYPE nucl) {
    int nuclValue = 4;
    switch (nucl) {
        case Params::A: {
            nuclValue = 0;
            break;
        }
        case Params::C: {
            nuclValue = 1;
            break;
        }
        case Params::G: {
            nuclValue = 2;
            break;
        }
        case Params::T: {
            nuclValue = 3;
            break;
        }
        case Params::N: {
            nuclValue = 15;
            break;
        }
        default: {
            cerr << "in getNuklNumber nucl = " << nucl << " - no such option, crashing..." << endl;
            cerr << char(nucl) << endl;
            exit(1);
        }

    }

    return nuclValue;

}

int Params::getNukl(char a) {
    if (a == 'A') return A;
    else if (a == 'C') return C;
    else if (a == 'G') return G;
    else if (a == 'T') return T;
    else return N;
}

string Params::getNuklAsString(int nukl) {
    string s = "";
    if (nukl == Params::A) s += 'A';
    else if (nukl == Params::C) s += 'C';
    else if (nukl == Params::G) s += 'G';
    else if (nukl == Params::T) s += 'T';
    else if (nukl == N) s += 'N';
    else {
        cerr << "returning \"\" as astring for nukl " << nukl << " in getNuklAsString" << endl;
        exit(1);
        return s;
    }

    return s;
}


void Params::initializeParams(int argc, char **argv) {
    long long defaultMaxHash = (100'000'000ll * 100'000'000ll + 1); // 10^16 + 1




    string file1 = "file1";
    string file2 = "file2";
    string alg = "alg";
    string remove_pref_reads = "remove_pref_reads";
    string rpr = "rpr";
    string li_kmer_intervals = "li_kmer_intervals";
    string li_kmer_length = "li_kmer_length";
    string redirect_cerr = "redirect_cerr";
    string max_offset_considered_for_alignment = "max_offset_considered_for_alignment";
    string mocfo = "mocfo";
    string max_offset_parallel_paths = "max_offset_parallel_paths";
    string mopp = "mopp";
    string max_offset_dangling_branches = "max_offset_dangling_branches";
    string modb = "modb";
    string contig_min_output_length = "contig_min_output_length";
    string max_error_rate_for_lcs = "max_error_rate_for_lcs";
    string minimal_overlap_rate_for_lcs = " minimal_overlap_rate_for_lcs";
    string minimal_overlap_for_lcs_low_error = "minimal_overlap_for_lcs_low_error";
    string min_overlap_rate = "min_overlap_rate";
    string mor = "mor";
    string min_overlap_area = "min_overlap_area";
    string moa = "moa";
    string use_lcs_low_error_filter = "use_lcs_low_error_filter";
    string use_acler_instead_of_aclcs = "use_acler_instead_of_aclcs";
    string acler = "acler";
    string offset_max_for_bitmap_alignment = "offset_max_for_bitmap_alignment";
    string serialize_graph_before_simplifier = "serialize_graph_before_simplifier";
    string serialize_graph_after_simplifier = "serialize_graph_after_simplifier";
    string serialize = "serialize";
    string deserialize_graph = "deserialize_graph";
    string threads = "threads";
    string write_statistics = "write_statistics";
    string new_reads_per_contig_percentage = "new_reads_per_contig_percentage";
    string nrpcp = "nrpcp";
    string alignment_controller_same_ends_length = "alignment_controller_same_ends_length";
    string acsel = "acsel";
    string rsoemo = "rsoemo";
    string rsoentr = "rsoentr";
    string read_end_trim = "read_end_trim";
    string read_end_trim_left = "read_end_trim_left";
    string retl = "retl";
    string read_end_trim_right = "read_end_trim_right";
    string retr = "retr";
    string remove_reads_with_n = "remove_reads_with_n";
    string min_offset_for_alignment = "min_offset_for_alignment";
    string rna = "rna";
    string scale = "scale";
    string output = "output";
    string tpn = "tpn";
    string tpv = "tpv";
    string error_rate = "error_rate";
    string er = "er";
    string correct_reads = "correct_reads";

    string mfup = "mfup"; // Most Frequently Used Parameter

//    int opt;
    string /*inStreamFilePath1,*/ fileName, fileNameNoExt;
    auto it = string::npos;
    int option_index = 0;

    static struct option long_options[] = {
            {file1.c_str(),         required_argument, 0, 0},
            {file2.c_str(),         required_argument, 0, 0},
            {threads.c_str(),       required_argument, 0, 0},
            {output.c_str(),        required_argument, 0, 0},
            {error_rate.c_str(),    required_argument, 0, 0},

//        {alg.c_str(),                                                   required_argument,          0,  0 },
//        {li_kmer_intervals.c_str(),                                     required_argument,          0,  0 },
//        {li_kmer_length.c_str(),                                        required_argument,          0,  0 },
            {redirect_cerr.c_str(), required_argument, 0, 0},
//        {max_offset_considered_for_alignment.c_str(),                   required_argument,          0,  0 },
//        {max_offset_parallel_paths.c_str(),                             required_argument,          0,  0 },
//        {max_offset_dangling_branches.c_str(),                          required_argument,          0,  0 },
//        {max_error_rate_for_lcs.c_str(),                                required_argument,          0,  0 },
//        {minimal_overlap_rate_for_lcs.c_str(),                          required_argument,          0,  0 },
//        {contig_min_output_length.c_str(),                              required_argument,          0,  0 },
//        {minimal_overlap_for_lcs_low_error.c_str(),                     required_argument,          0,  0 },
//        {use_lcs_low_error_filter.c_str(),                              required_argument,          0,  0 },
//        {use_acler_instead_of_aclcs.c_str(),                            required_argument,          0,  0 },
//        {offset_max_for_bitmap_alignment.c_str(),                       required_argument,          0,  0 },
//        {serialize_graph_before_simplifier.c_str(),                     required_argument,          0,  0 },
//        {serialize_graph_after_simplifier.c_str(),                      required_argument,          0,  0 },
            {deserialize_graph.c_str(), required_argument, 0, 0},
//        {min_overlap_rate.c_str(),                                      required_argument,          0,  0 },
//        {write_statistics.c_str(),                                      required_argument,          0,  0 },
//        {mocfo.c_str(),                                                 required_argument,          0,  0 },
//        {mor.c_str(),                                                   required_argument,          0,  0 },
//        {acler.c_str(),                                                 required_argument,          0,  0 },
//        {mopp.c_str(),                                                  required_argument,          0,  0 },
//        {modb.c_str(),                                                  required_argument,          0,  0 },
//        {mfup.c_str(),                                                  required_argument,          0,  0 },
//        {nrpcp.c_str(),                                                 required_argument,          0,  0 },
//        {new_reads_per_contig_percentage.c_str(),                       required_argument,          0,  0 },
//        {alignment_controller_same_ends_length.c_str(),                 required_argument,          0,  0 },
//        {acsel.c_str(),                                                 required_argument,          0,  0 },
//        {min_overlap_area.c_str(),                                      required_argument,          0,  0 },
//        {rsoemo.c_str(),                                                required_argument,          0,  0 },
//        {moa.c_str(),                                                   required_argument,          0,  0 },
//        {read_end_trim.c_str(),                                         required_argument,          0,  0 },
//        {read_end_trim_left.c_str(),                                    required_argument,          0,  0 },
//        {read_end_trim_right.c_str(),                                   required_argument,          0,  0 },
//        {retl.c_str(),                                                  required_argument,          0,  0 },
//        {retr.c_str(),                                                  required_argument,          0,  0 },
//        {rsoentr.c_str(),                                               required_argument,          0,  0 },
            {remove_reads_with_n.c_str(), required_argument, 0, 0},
//        {min_offset_for_alignment.c_str(),                              required_argument,          0,  0 },
//        {rpr.c_str(),                                                   required_argument,          0,  0 },
//        {remove_pref_reads.c_str(),                                     required_argument,          0,  0 },
            {serialize.c_str(),           required_argument, 0, 0},
            {rna.c_str(),                 required_argument, 0, 0},
            {scale.c_str(),               required_argument, 0, 0},
//        {tpn.c_str(),                                                required_argument,          0,  0 },
//        {tpv.c_str(),                                                required_argument,          0,  0 },
//        {er.c_str(),                                                required_argument,          0,  0 },
//        {correct_reads.c_str(),                                                required_argument,          0,  0 },
            {0, 0,                                           0, 0}
    };


    while (1) {

        int this_option_optind = optind ? optind : 1;
        int option_index = 0;
        int c;
        string option, option_name;

        c = getopt_long(argc, argv, "l:",
                        long_options, &option_index);
        if (c == -1) break;
        switch (c) {
            case 0:
                option = string(optarg);
                option_name = string(long_options[option_index].name);
//                cerr << "\toption = " << option << "    option_name = " << option_name << endl;

                if (option_name == file1) {
                    inStreamFilePath1 = string(optarg);

                    cerr << "first inStreamFilePath1: " << inStreamFilePath1 << endl;
//                    inStream.open(inStreamFilePath1);
//                    cin.rdbuf(inStream.rdbuf());

//                    fileName;
                    it = inStreamFilePath1.rfind('/');
                    if (it == string::npos) fileName = inStreamFilePath1;
                    else fileName = inStreamFilePath1.substr(it + 1);


                    fileExtension = "";
                    it = fileName.rfind('.');
                    if (it == string::npos) fileExtension = "";
                    else fileExtension = fileName.substr(it + 1);


                    if (fileExtension == "fasta") INPUT_FILE_TYPE = FASTA;
                    else if (fileExtension == "pfasta") INPUT_FILE_TYPE = PFASTA;
                    else if (fileExtension == "fastq") INPUT_FILE_TYPE = FASTQ;
                    else INPUT_FILE_TYPE = MY_INPUT;

                    fileNameNoExt = "";
                    if (fileExtension == "") fileNameNoExt = fileName;
                    else fileNameNoExt = fileName.substr(0, fileName.size() - fileExtension.size() - 1);


//                    if( fileNameNoExt != "" ) fileNameNoExt.erase( fileNameNoExt.size()-1,1 );
                    TEST_NAME = "ALGA_" + fileNameNoExt;

                }
                if (option_name == error_rate || option_name == er) {
                    double rate = stod(string(optarg));
                    /* if( rate <= 0.01 ) ALGORITHM_IN_USE = PREF_SUF_GRAPH_CREATION;
                     else{
                         ALGORITHM_IN_USE = PAIRWISE_KMER_BRANCH_GRAPH_CREATION;
                         if( rate < 0.05 ){
                             MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = 95;
                         }else{
                             MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = (1-rate) * 100;
                         }
                     }*/
                    ERROR_RATE = 100 * rate;
                    if (rate <= 0.01) USE_GRAPH_CREATOR_SUPPLEMENT = 0;
                    else USE_GRAPH_CREATOR_SUPPLEMENT = 1;
                }
                if (option_name == output) {
                    outStreamFileName = string(optarg);
                }
                if (option_name == correct_reads) {
                    CORRECT_READS = stoi(string(optarg));
                }
                if (option_name == tpn) {
                    TPN = string(optarg);
                }
                if (option_name == tpv) {
                    TPV = string(optarg);
                }
                if (option_name == file2) {
                    inStreamFilePath2 = string(optarg);
                }
                if (option_name == rpr || option_name == remove_pref_reads) {
                    string opt = string(optarg);
                    if (opt == "all" || opt == "a") REMOVE_PREF_READS_TYPE = PREF_READS_ALL_PREFIX_READS;
                    else if (opt == "duplicates" || opt == "d") REMOVE_PREF_READS_TYPE = PREF_READS_ONLY_DUPLICATES;
                    else if (opt == "none" || opt == "n") REMOVE_PREF_READS_TYPE = PREF_READS_NONE;
                }
                if (option_name == scale) {
                    SCALE = stof(string(optarg));
                }
                if (option_name == alg) {
                    if (string(optarg) == "PKB") ALGORITHM_IN_USE = PAIRWISE_KMER_BRANCH_GRAPH_CREATION;
                    else if (string(optarg) == "PS") ALGORITHM_IN_USE = PREF_SUF_GRAPH_CREATION;
                }
                if (option_name == rna) {
                    RNA = stoi(string(optarg));
                }
                if (option_name == serialize) {
                    SERIALIZE_GRAPH_AFTER_SIMPLIFIER = SERIALIZE_GRAPH_BEFORE_SIMPLIFIER = DESERIALIZE_GRAPH = stoi(
                            string(optarg));
                }
                if (option_name == li_kmer_intervals) {
                    LI_KMER_INTERVALS = stoi(string(optarg));
                }
                if (option_name == li_kmer_length) {
                    LI_KMER_LENGTH = stoi(string(optarg));
                }
                if (option_name == redirect_cerr) {
                    REDIRECT_CERR = stoi(string(optarg));
                }
                if (option_name == max_offset_considered_for_alignment || option_name == mocfo) {
                    MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT = stoi(string(optarg));
                }
                if (option_name == max_offset_parallel_paths || option_name == mopp) {
                    MAX_OFFSET_PARALLEL_PATHS = stoi(string(optarg));
                }
                if (option_name == max_offset_dangling_branches || option_name == modb) {
                    MAX_OFFSET_DANGLING_BRANCHES = stoi(string(optarg));
                }
                if (option_name == max_error_rate_for_lcs) {
                    MAX_ERROR_RATE_FOR_LCS = stoi(string(optarg));
                }
                if (option_name == minimal_overlap_rate_for_lcs) {
                    MINIMAL_OVERLAP_RATE_FOR_LCS = stoi(string(optarg));
                }
                if (option_name == contig_min_output_length) {
                    CONTIG_MIN_OUTPUT_LENGTH = stoi(string(optarg));
                }
                if (option_name == minimal_overlap_for_lcs_low_error) {
                    MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = stoi(string(optarg));
                }
                if (option_name == use_lcs_low_error_filter) {
                    USE_LCS_LOW_ERROR_FILTER = stoi(string(optarg));
                }
                if (option_name == use_acler_instead_of_aclcs || option_name == acler) {
                    USE_ACLER_INSTEAD_OF_ACLCS = stoi(string(optarg));
                }
                if (option_name == serialize_graph_before_simplifier) {
                    SERIALIZE_GRAPH_BEFORE_SIMPLIFIER = stoi(string(optarg));
                }
                if (option_name == serialize_graph_after_simplifier) {
                    SERIALIZE_GRAPH_AFTER_SIMPLIFIER = stoi(string(optarg));
                }
                if (option_name == mfup) {
                    MOST_FREQUENTLY_USED_PARAMETER = stoi(string(optarg));
                    KMER_LENGTH_BUCKET = stoi(string(optarg));
                    LI_KMER_LENGTH = stoi(string(optarg));
                    MIN_OVERLAP_PREF_SUF = stoi(string(optarg));

//                    cerr << "MIN_OVERLAP_PREF_SUF = " << MIN_OVERLAP_PREF_SUF << endl;

                    MIN_OVERLAP_AREA = stoi(string(optarg));
                }
                if (option_name == min_overlap_rate || option_name == mor) {
                    MIN_OVERLAP_RATE = stoi(string(optarg));
                    MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = (100 + stoi(string(optarg))) >> 1;
                    MINIMAL_OVERLAP_RATE_FOR_LCS = stoi(string(optarg));
                }
                if (option_name == deserialize_graph) {
                    DESERIALIZE_GRAPH = stoi(string(optarg));
                }
                if (option_name == write_statistics) {
                    WRITE_STATISTICS = stoi(string(optarg));
                }
                if (option_name == threads) {
                    THREADS = stoi(string(optarg));
                }
                if (option_name == new_reads_per_contig_percentage || option_name == nrpcp) {
                    NEW_READS_PER_CONTIG_PERCENTAGE = stoi(string(optarg));
                }
                if (option_name == alignment_controller_same_ends_length || option_name == acsel) {
                    ALIGNMENT_CONTROLLER_SAME_ENDS_LENGTH = stoi(string(optarg));
                }
                if (option_name == min_overlap_area || option_name == moa) {
                    MIN_OVERLAP_AREA = stoi(string(optarg));
                }
                if (option_name == rsoemo) {
                    REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP = stoi(string(optarg));
                }
                if (option_name == rsoentr) {
                    REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN = stoi(string(optarg));
                }
                if (option_name == read_end_trim) {
                    READ_END_TRIM_LEFT = READ_END_TRIM_RIGHT = stoi(string(optarg));
                }
                if (option_name == read_end_trim_left || option_name == retl) {
                    READ_END_TRIM_LEFT = stoi(string(optarg));
                }
                if (option_name == read_end_trim_right || option_name == retr) {
                    READ_END_TRIM_RIGHT = stoi(string(optarg));
                }
                if (option_name == remove_reads_with_n) {
                    REMOVE_READS_WITH_N = stoi(string(optarg));
                }

                break;
            case 'l':
                MOST_FREQUENTLY_USED_PARAMETER = stoi(string(optarg));

                KMER_LENGTH_BUCKET = stoi(string(optarg));
                LI_KMER_LENGTH = stoi(string(optarg));
                MIN_OVERLAP_PREF_SUF = stoi(string(optarg));

                MIN_OVERLAP_AREA = stoi(string(optarg));

                break;
            case '?':
                break;
            default:
                printf("?? getopt returned character code 0%o ??\n", c);
        }
    }


    std::time_t end_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    string date = string(std::ctime(&end_time));
    date = MyUtils::replaceAll(date, " ", "_");
    date = MyUtils::replaceAll(date, ":", "_");
    date = MyUtils::replaceAll(date, "\n", "");
    RUN_DATE = date;


    /*outStreamFileName = "";
    if( USE_DATE_IN_STREAM_FILE_NAMES ){
        outStreamFileName += "_" + date;
    }

    if( ALGORITHM_IN_USE == PREF_SUF_GRAPH_CREATION ){
        TEST_NAME += "_alg" + to_string( ALGORITHM_IN_USE ) + "_l" + ( MIN_OVERLAP_PREF_SUF == -1 ? "def" : to_string( MIN_OVERLAP_PREF_SUF ) );
        TEST_NAME += "_rpr" + to_string(REMOVE_PREF_READS_TYPE);
            TEST_NAME += "_rsoe" + ( REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP == -1 ? "def" : to_string( REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP ) )
                    + "-" + to_string( Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN );

    }
    else{
        TEST_NAME += "_alg" + to_string( ALGORITHM_IN_USE );
        TEST_NAME += "_rpr" + to_string(REMOVE_PREF_READS_TYPE);
        TEST_NAME += "_LI-l" + to_string( LI_KMER_LENGTH ) + "-int" +to_string( LI_KMER_INTERVALS );

        TEST_NAME += "_rsoe" + ( REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP == -1 ? "def" : to_string( REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP ) )
                     + "-" + to_string( Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN );

        if( USE_ACLER_INSTEAD_OF_ACLCS ) TEST_NAME += "_acler";

        TEST_NAME += "_mocfo" + to_string(MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT);
        TEST_NAME += "_mor" + to_string( MIN_OVERLAP_RATE );
        TEST_NAME += "_moa" + to_string( MIN_OVERLAP_AREA );

        TEST_NAME += "_acsel" + to_string(ALIGNMENT_CONTROLLER_SAME_ENDS_LENGTH);
    }

//    if( ADD_COMP_REV_READS ) TEST_NAME += "_acrr";
//    if( ADD_PAIRED_READS ) TEST_NAME += "_apr";

    TEST_NAME += "_ret" + to_string(READ_END_TRIM_LEFT) + "-" + to_string(READ_END_TRIM_RIGHT); // + "_mofa" + to_string(MIN_OFFSET_FOR_ALIGNMENT);

    outStreamFileName += TEST_NAME;*/


    TEST_NAME += "_scale" + to_string((int) (100 *
                                             SCALE));// + "_mfup" + to_string(MOST_FREQUENTLY_USED_PARAMETER) + "_rsoemo" + to_string(REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP);

    TEST_NAME += (REMOVE_READS_WITH_N ? "_noN" : "_randN");

    cerr << "inStreamFilePath1 = " << inStreamFilePath1 << endl;
    cerr << "fileName = " << fileName << endl;
    cerr << "fileExtension = " << fileExtension << endl;
    cerr << "fileNameNoExt = " << fileNameNoExt << endl;
    cerr << "TEST_NAME = " << TEST_NAME << endl;


    { // 'master' version
        if (inStreamFilePath1 == "") {
            cerr << endl << "ERROR - PLEASE PROVIDE THE INPUT FILE using --file1 option!" << endl;
            cout << endl << "ERROR - PLEASE PROVIDE THE INPUT FILE using --file1 option!" << endl;
            exit(1);
        }

        if (outStreamFileName == "") {
            cerr << endl << "ERROR - PLEASE PROVIDE THE OUTPUT FILE NAME!" << endl;
            cout << endl << "ERROR - PLEASE PROVIDE THE OUTPUT FILE NAME!" << endl;
            exit(1);
        }
        {
            string logName = outStreamFileName;
            auto it = logName.find(".fasta");
            if (it == string::npos) it = logName.find(".fastq");
            if (it != string::npos) logName.erase(it);
            errStream.open(logName + ".log");
//            cerr.rdbuf(errStream.rdbuf());
        }
    }


    if (REDIRECT_CERR) {
        string tpn = "_tpn-" + TPN + "_tpv-" + TPV;
        errStream.open((USE_DATE_IN_STREAM_FILE_NAMES ? date : "") +
                       (outStreamFileName != "" ? outStreamFileName : TEST_NAME) + (TPN == "" ? "" : tpn) + ".log");
        cerr.rdbuf(errStream.rdbuf());
    }

}

void Params::createOutputFileName() {

    string tpn = "_tpn-" + TPN + "_tpv-" + TPV;

    if (outStreamFileName != "") {
        if (TPN != "") outStreamFileName += tpn;

        /*if( REDIRECT_CERR ){
            errStream.close();
            errStream.open( outStreamFileName + ".log" );
            cerr.rdbuf( errStream.rdbuf() );
        }*/

        cerr << "outStreamFileName: " << outStreamFileName << endl;
        return;
    }
    outStreamFileName = "";


    std::time_t end_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    string date = string(std::ctime(&end_time));
    date = MyUtils::replaceAll(date, " ", "_");
    date = MyUtils::replaceAll(date, ":", "_");
    date = MyUtils::replaceAll(date, "\n", "");
    RUN_DATE = date;

    if (USE_DATE_IN_STREAM_FILE_NAMES) {
        outStreamFileName += "_" + date;
    }

    if (TPN != "") TEST_NAME += tpn;

//    TEST_NAME += "_scale" + to_string((int)(100*SCALE));
    if (ALGORITHM_IN_USE == PREF_SUF_GRAPH_CREATION) {
        TEST_NAME += "_alg" + to_string(ALGORITHM_IN_USE) + "_l" +
                     (MIN_OVERLAP_PREF_SUF == -1 ? "def" : to_string(MIN_OVERLAP_PREF_SUF));
        TEST_NAME += "_rpr" + to_string(REMOVE_PREF_READS_TYPE);
        TEST_NAME += "_rsoe" + (REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP == -1 ? "def" : to_string(
                REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP))
                     + "-" + to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN);

    } else {
        TEST_NAME += "_alg" + to_string(ALGORITHM_IN_USE);
        TEST_NAME += "_rpr" + to_string(REMOVE_PREF_READS_TYPE);
        TEST_NAME += "_LI-l" + to_string(LI_KMER_LENGTH) + "-int" + to_string(LI_KMER_INTERVALS);

        TEST_NAME += "_rsoe" + (REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP == -1 ? "def" : to_string(
                REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP))
                     + "-" + to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN);

        if (USE_ACLER_INSTEAD_OF_ACLCS) TEST_NAME += "_acler";

        TEST_NAME += "_mocfo" + to_string(MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT);
        TEST_NAME += "_mor" + to_string(MIN_OVERLAP_RATE);
        TEST_NAME += "_moa" + to_string(MIN_OVERLAP_AREA);

        TEST_NAME += "_acsel" + to_string(ALIGNMENT_CONTROLLER_SAME_ENDS_LENGTH);
    }


//    if( ADD_COMP_REV_READS ) TEST_NAME += "_acrr";
//    if( ADD_PAIRED_READS ) TEST_NAME += "_apr";

    TEST_NAME += "_ret" + to_string(READ_END_TRIM_LEFT) + "-" +
                 to_string(READ_END_TRIM_RIGHT); // + "_mofa" + to_string(MIN_OFFSET_FOR_ALIGNMENT);


    outStreamFileName += TEST_NAME;


    cerr << "TEST_NAME = " << TEST_NAME << endl;

    /*if( REDIRECT_CERR ){
        errStream.open( outStreamFileName + ".log" );
        cerr.rdbuf( errStream.rdbuf() );
    }*/
}


int Params::MOST_FREQUENTLY_USED_PARAMETER = -1; // instead of     KMER_LENGTH_BUCKET,     LI_KMER_LENGTH     and     MIN_OVERLAP_PREF_SUF
float Params::SCALE = 0.55;


int Params::KMER_LENGTH_BUCKET = MOST_FREQUENTLY_USED_PARAMETER;


int Params::MAX_OFFSET_CONSIDERED_FOR_ALIGNMENT = 70;


int Params::MAX_OFFSET_PARALLEL_PATHS = 250;
int Params::MAX_OFFSET_DANGLING_BRANCHES = 250;


int Params::WRITE_STATISTICS = 0;
int Params::ADD_COMP_REV_READS = 1;
int Params::ADD_PAIRED_READS = 1;


int Params::MIN_OVERLAP_RATE = 95;

int Params::MINIMAL_OVERLAP_RATE_FOR_LCS = MIN_OVERLAP_RATE;
int Params::MAX_ERROR_RATE_FOR_LCS = 2;

int Params::MINIMAL_OVERLAP_FOR_LCS_LOW_ERROR = (100 + MIN_OVERLAP_RATE) >> 1;
int Params::USE_LCS_LOW_ERROR_FILTER = 1;
int Params::USE_ACLER_INSTEAD_OF_ACLCS = 1;

int Params::LI_KMER_LENGTH = MOST_FREQUENTLY_USED_PARAMETER;
int Params::LI_KMER_INTERVALS = 3;

int Params::MIN_OVERLAP_PREF_SUF = MOST_FREQUENTLY_USED_PARAMETER;
int Params::MIN_OFFSET_FOR_ALIGNMENT = 0;
int Params::MIN_OVERLAP_AREA = -1;


int Params::PAIRWISE_KMER_BRANCH_GRAPH_CREATION = 3;
int Params::PREF_SUF_GRAPH_CREATION = 5;


int Params::ALGORITHM_IN_USE = PREF_SUF_GRAPH_CREATION;


//Params::KMER_HASH_TYPE Params::MAX_HASH_CONSIDERED = ((100000000ll * 100000000ll + 1));
Params::KMER_HASH_TYPE Params::MAX_HASH_CONSIDERED = ((1'000'000'000ll * 1'000'000'000ll + 3)); // 10^18 + 3


int Params::TRAVERSE_TYPE = /*TRAVERSE_GREEDY;*/    TRAVERSE_SHALLOW;


int Params::REMOVE_SHORT_PARALLEL_PATHS_MST_FAST = 1;

int Params::READ_END_TRIM_LEFT = 3;
int Params::READ_END_TRIM_RIGHT = 3;

int Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP = -1;
int Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN = 3;


int Params::CONTIG_MIN_OUTPUT_LENGTH = 200;
int Params::CONTIG_CREATOR_SHORT_CYCLE_LENGTH = 150;


int Params::INPUT_FILE_TYPE = MY_INPUT;
int Params::REMOVE_READS_WITH_N = 1;

int Params::CORRECT_READS = 0;


bool Params::REDIRECT_CERR = false;

int Params::SERIALIZE_GRAPH_BEFORE_SIMPLIFIER = 0;
int Params::SERIALIZE_GRAPH_AFTER_SIMPLIFIER = 0;
int Params::DESERIALIZE_GRAPH = 0;

string Params::TPN = "";
string Params::TPV = "";

int Params::NEW_READS_PER_CONTIG_PERCENTAGE = 95;
int Params::ALIGNMENT_CONTROLLER_SAME_ENDS_LENGTH = 3;

int Params::REMOVE_PREF_READS_TYPE = PREF_READS_ALL_PREFIX_READS;

int Params::USE_GRAPH_CREATOR_SUPPLEMENT = 0;

int Params::LI_PRIORITIES_TO_CONSIDER = 4;
int Params::THREADS = 6;

string Params::TEST_NAME = "ALGA_NO_NAME_TESTING";
string Params::fileExtension = "";
string Params::outStreamFileName = "";
string Params::inStreamFilePath1 = "";
string Params::inStreamFilePath2 = "";
int Params::USE_DATE_IN_STREAM_FILE_NAMES = 0;
string Params::RUN_DATE = "NO_DATE";

int Params::RNA = 0;
int ERROR_RATE = 0;

ifstream Params::inStream;
ofstream Params::outStream;
ofstream Params::errStream;