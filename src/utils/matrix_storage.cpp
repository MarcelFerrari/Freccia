#include <Eigen/Dense>
#include "Freccia/utils/matrix_storage.hpp"

// Function to convert a matrix to the banded compact storage format
Eigen::MatrixXd Freccia::Banded::make_compressed(const Eigen::Ref<const Eigen::MatrixXd> & A, const unsigned int f){
    unsigned int n = A.cols();

    Eigen::MatrixXd A_band = Eigen::MatrixXd::Zero(f, n);
   
    // Copy matrix in banded format
    for(unsigned j = 0; j < n; ++j){
        for(unsigned i = j; i < std::min(n, j+f); ++i){
            A_band(i-j, j) = A(i,j);
        }
    }

    return A_band;
};

// Function to convert a banded matrix to the dense format
Eigen::MatrixXd Freccia::Banded::make_dense(const Eigen::Ref<const Eigen::MatrixXd> & A_band, const unsigned int f){
    unsigned int n = A_band.cols();  // assuming the banded matrix has the same number of columns as the original matrix

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);  // initialize a zero matrix of size n x n

    // Copy the banded matrix back to the dense format
    for(unsigned j = 0; j < n; ++j){
        for(unsigned i = j; i < std::min(n, j+f); ++i){
            A(i, j) = A_band(i-j, j);
        }
    }

    return A;
}
