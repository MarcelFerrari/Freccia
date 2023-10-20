#ifndef FRECCIA_BANDED_EIGENSOLVER_H
#define FRECCIA_BANDED_EIGENSOLVER_H
#include <Eigen/Dense>
#include <iostream>
#include <vector>

// Wrapper class for LAPACK dsbevd routine
namespace Freccia::Banded {

class BandedEigenSolver {
public:
    BandedEigenSolver(const Eigen::Ref<const Eigen::MatrixXd> & A, const unsigned int f_in, char uplo_in = '\0')
    :f(f_in), uplo(uplo_in), ew(A.cols()), ev(A.cols(), A.cols())
    {
        // Convert matrix to the banded storage format if necessary
        if (uplo != '\0') {
            A_band = A;
        } else {
            make_banded(A);  // Assuming this modifies 'A' in-place or updates some member of the class
        }

        // Perform eigendecomposition
        dsbevd_wrapper();
    }

    // Getter functions
    const Eigen::VectorXd& eigenvalues() const { return ew; }
    const Eigen::MatrixXd& eigenvectors() const { return ev; }

private:
    unsigned int f;
    char uplo; // Assuming matrix is stored in lower banded format
    Eigen::MatrixXd A_band;  // The matrix in banded format
    Eigen::VectorXd ew;   // Eigenwerte
    Eigen::MatrixXd ev;  // Eigenvektoren

    void make_banded(const Eigen::Ref<const Eigen::MatrixXd> & A);
    void dsbevd_wrapper();
};
}
#endif