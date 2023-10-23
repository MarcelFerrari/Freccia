#ifndef FRECCIA_MATRIX_STORAGE_H
#define FRECCIA_MATRIX_STORAGE_H

#include <Eigen/Dense>

namespace Freccia::Banded {

// Function to convert a matrix to the banded compact storage format
Eigen::MatrixXd make_compressed(const Eigen::Ref<const Eigen::MatrixXd> & A, const unsigned int f);

}

#endif