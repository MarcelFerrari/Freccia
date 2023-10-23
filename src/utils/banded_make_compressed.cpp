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