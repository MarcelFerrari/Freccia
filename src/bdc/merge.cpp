// Standard libraries
#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <tuple>
#include <limits>

// Include eigen
#include <Eigen/Dense>

// Custom headers
#include "Freccia/bdc/BDCEigenSolver.hpp"
#include "Freccia/dprk/DPRKEigenSolver.hpp"

// Convenience aliases
using SVDTuple = std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd>;
using EigenTuple = std::tuple<unsigned int, unsigned int>;

EigenTuple Freccia::Banded::BDCEigenSolver::mergeBlocks(EigenTuple& A1, EigenTuple& A2){
    // Get start and size of blocks
    unsigned int start1 = std::get<0>(A1);
    unsigned int size1 = std::get<1>(A1);
    unsigned int start2 = std::get<0>(A2);
    unsigned int size2 = std::get<1>(A2);

    // Get references to blocks
    Eigen::Ref<Eigen::VectorXd> D_ = D.segment(start1*f, (size1 + size2)*f);
    Eigen::Ref<Eigen::MatrixXd> Q_ = Q.block(start1*f, start1*f, (size1 + size2)*f, (size1 + size2)*f);

    // Get required SVD matrices
    const Eigen::MatrixXd& U = std::get<0>(svdBlocks[start1 + size1 - 1]);
    const Eigen::VectorXd& S = std::get<1>(svdBlocks[start1 + size1 - 1]);
    const Eigen::MatrixXd& V = std::get<2>(svdBlocks[start1 + size1 - 1]);

    // Compute W
    // W is N*q where q is the number of columns of V*S or U*S
    unsigned int q = V.cols(); // V.cols() == U.cols()
    
    // Allocate memory for W
    Eigen::MatrixXd W = Eigen::MatrixXd::Zero((size1 + size2)*f, q);
    Eigen::VectorXd sqrtS = S.array().sqrt();
    
    W.middleRows(size1*f - q, q) = V*sqrtS.asDiagonal();
    W.middleRows(size1*f, q) = U*sqrtS.asDiagonal();

    // Solve DPRK problem
    Freccia::DPRK::DPRKEigenSolver solver(D_, W, 1., Q_); // Eigensystem update problem
    D_ = solver.eigenvalues();
    Q_ = solver.eigenvectors();

    return std::make_tuple(start1, size1 + size2); // Encode result
}