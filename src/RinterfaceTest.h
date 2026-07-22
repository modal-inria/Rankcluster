#ifndef RUNTEST_H_
#define RUNTEST_H_

#include <RcppEigen.h>
#include <vector>
#include <set>
#include "test.h"

std::vector<std::vector<double>> convertToVVd(Rcpp::NumericMatrix data);
std::vector<std::vector<int>> convertToVVi(Rcpp::NumericMatrix data);
std::vector<RankStruct> downUniVariateRank(Rcpp::NumericMatrix XR);

#endif /* RUNTEST_H_ */
