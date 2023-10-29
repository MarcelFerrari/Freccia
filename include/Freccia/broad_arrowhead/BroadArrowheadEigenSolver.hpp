#ifndef FRECCIA_ARROWHEAD_BROADARROWHEADEIGENSOLVER_H
#define FRECCIA_ARROWHEAD_BROADARROWHEADEIGENSOLVER_H

// Eigen
#include <Eigen/Dense>
#include <string>

namespace Freccia::Arrowhead {

    class BroadArrowheadEigenSolver {
        public:
        // Constructor for compact matrix storage
        // B is assumed to be in compact stoarage format
        BroadArrowheadEigenSolver(const Eigen::Ref<const Eigen::MatrixXd>& B, const Eigen::Ref<const Eigen::MatrixXd>& W, unsigned int f_in, unsigned int g_in, const std::string & method_in="banded", char uplo_in = 'L')
        : f(f_in), g(g_in), uplo(uplo_in), l(B.cols()), n(l+g), method(method_in) {
            eigh(B, W);
        };

        // Getters
        const Eigen::VectorXd& eigenvalues() const { return ew; };
        
        const Eigen::MatrixXd& eigenvectors() const { return ev; };
        
        private:
        // Member variables
        const unsigned int l; // Size of the banded block
        const unsigned int f; /* If the banded eigensolver is used, this is the bandwidth of the shaft.
                               * If the BDC algorithm is used, this is the size of the diagonal and off-diagonal blocks. */
        const char uplo;      // Storage format of the shaft
        const unsigned int g; // Width of the arrowhead
        const unsigned int n; // Size of the matrix
        const std::string& method; // Method used to solve the eigenvalue problem

        // Result
        Eigen::VectorXd ew; // Eigenwerte
        Eigen::MatrixXd ev; // Eigenvektoren

        // Member functions
        void eigh(const Eigen::Ref<const Eigen::MatrixXd>& B, const Eigen::Ref<const Eigen::MatrixXd>& W); // Main driver
    };

} // End namespace Freccia::DPRK
#endif // FRECCIA_ARROWHEAD_BROADARROWHEADEIGENSOLVER_H