//
// Created by sylwester on 12/20/18.
//

#include "StatisticsGenerators/StatisticsGeneratorBigData.h"

string StatisticsGeneratorBigData::SAME_KMER_SIZES = "SAME_KMER_SIZES";
string StatisticsGeneratorBigData::SWAT_ALG_BRANCHES = "BRANCHES IN SWAT ALG";

string StatisticsGeneratorBigData::COMPARISONS_IN_PAIRWISE_KMER_ALGORITHM = "READ COMPARISONS_IN_PAIRWISE_KMER_ALGORITHM";
string StatisticsGeneratorBigData::POSITIVE_COMPARISONS_IN_PAIRWISE_KMER_ALGORITHM = "READ COMPARISONS_IN_PAIRWISE_KMER_ALGORITHM - POSITIVE";
string StatisticsGeneratorBigData::COMPARISONS_IN_SWAT_BRANCH_ALGORITHM = "READ COMPARISONS_IN_SWAT_BRANCH_ALGORITHM";
map<string, StatisticsGeneratorBigData::Stat> StatisticsGeneratorBigData::data = map<string, StatisticsGeneratorBigData::Stat>();