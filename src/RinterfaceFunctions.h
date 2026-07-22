#ifndef RUNFUNCTIONS_H_
#define RUNFUNCTIONS_H_

#include <RcppEigen.h>
#include <vector>

std::vector<std::vector<std::vector<int>>> numMat2vvvInt(Rcpp::NumericMatrix XR, std::vector<int> const &m);

#endif /* RUNFUNCTIONS_H_ */
