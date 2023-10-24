#ifndef FRECCIA_ARROWHEAD_BROADARROWHEADEIGENSOLVER_H
#define FRECCIA_ARROWHEAD_BROADARROWHEADEIGENSOLVER_H

// Eigen
#include <Eigen/Dense>

namespace Freccia::Arrowhead {

    class BroadArrowheadEigenSolver {
        public:
        // Constructor for compact matrix storage
        // B is assumed to be in compact stoarage format
        BroadArrowheadEigenSolver(const Eigen::Ref<const Eigen::MatrixXd>& B, const Eigen::Ref<const Eigen::MatrixXd>& W, char uplo_in)
        : f(B.rows()), g(W.cols()), uplo(uplo_in), n(B.cols() + W.cols()), l(B.cols())
        {
            eigh(B, W);
        };

        // Getters
        const Eigen::VectorXd& eigenvalues() const { return ew; };
        
        const Eigen::MatrixXd& eigenvectors() const { return ev; };
        
        private:
        // Member variables
        unsigned int l; // Size of the banded block
        unsigned int f; // Bandwidth of the shaft
        char uplo;      // Storage format of the shaft
        unsigned int g; // Bandwidth of the arrow head
        unsigned int n; // Size of the matrix

        // Result
        Eigen::VectorXd ew; // Eigenwerte
        Eigen::MatrixXd ev; // Eigenvektoren

        // Member functions
        void eigh(const Eigen::Ref<const Eigen::MatrixXd>& B, const Eigen::Ref<const Eigen::MatrixXd>& W); // Main driver
    };

} // End namespace Freccia::DPRK
#endif // FRECCIA_ARROWHEAD_BROADARROWHEADEIGENSOLVER_H