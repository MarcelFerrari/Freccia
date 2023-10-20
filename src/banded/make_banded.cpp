#include <Eigen/Dense>
#include "Freccia/banded/BandedEigenSolver.hpp"
#include <iostream>

// Function to convert a matrix to the banded compact storage format
void Freccia::Banded::BandedEigenSolver::make_banded(const Eigen::Ref<const Eigen::MatrixXd> & A){
    unsigned int n = A.cols();

    A_band = Eigen::MatrixXd::Zero(f, n);
   
    // Copy matrix in banded format
    for(unsigned j = 0; j < n; ++j){
        for(unsigned i = j; i < std::min(n, j+f); ++i){
            A_band(i-j, j) = A(i,j);
        }
    }
    
    // Set uplo to L storage
    uplo = 'L';
};