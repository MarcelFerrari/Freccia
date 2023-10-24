#include <Eigen/Dense>
#include <iostream>

// Freccia
#include "Freccia/broad_arrowhead/BroadArrowheadEigenSolver.hpp"
#include "Freccia/arrowhead/ArrowheadEigenSolver.hpp"
#include "Freccia/banded/BandedEigenSolver.hpp"

// Full eigenproblem
void Freccia::Arrowhead::BroadArrowheadEigenSolver::eigh(const Eigen::Ref<const Eigen::VectorXd> & B, const Eigen::Ref<const Eigen::MatrixXd> & W) {
    // Initialize result
    ew = Eigen::VectorXd::Zero(n);
    ev = Eigen::MatrixXd::Zero(n, n);

    // Diagonalize banded block
    { // Do not pollute namespace
        Freccia::Banded::BandedEigenSolver solver(B, f, uplo);
        ew.head(l) = solver.eigenvalues();
        ev.topLeftCorner(l, l) = solver.eigenvectors();
    }

    // Diagonalize using explicit form of recursion
    for(unsigned int k = 0; k < g; ++k){
        // Compute updated arrowhead matrix
        Eigen::VectorXd u = ev.topLeftCorner(l+k, l+k).transpose() * W.col(k).head(l+k);
        double alpha = W(l+k, k);
        
        // Diagonalize arrowhead matrix
        Freccia::Arrowhead::ArrowheadEigenSolver solver(ew.head(l+k), u, alpha);
        
        // Update eigenvalues and eigenvectors
        ew.head(l+k+1) = solver.eigenvalues();
        ev.topLeftCorner(l+k+1, l+k+1) *= solver.eigenvectors(); // Backpropagation of eigenvectors
    }

}
