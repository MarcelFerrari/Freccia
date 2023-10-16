#include <Eigen/Dense>
#include <iostream>
// Freccia
#include "Freccia/dprk/DPRKEigenSolver.hpp"
#include "Freccia/dpr1/DPR1EigenSolver.hpp"

// Full eigenproblem
void Freccia::DPRK::DPRKEigenSolver::DPRKeigh(const Eigen::Ref<const Eigen::VectorXd> & D, const Eigen::Ref<const Eigen::MatrixXd> & W, const double rho) {
    // Get dimensions of problem
    assert(D.size() == W.rows());
    unsigned int N = W.rows();
    unsigned int f = W.cols();


    // Solve first subproblem    
    { // Local scope to avoid polluting namespace
        Eigen::VectorXd z = W.col(0);
        
        Freccia::DPR1::DPR1EigenSolver eigh(D, z, rho);
        
        ew = eigh.eigenvalues();
        ev = eigh.eigenvectors();
    }

    // Compute DPRK update
    for(int i = 1; i < f; ++i){
        // Compute i-th DPR1 update
        Eigen::VectorXd z = ev.transpose()*W.col(i);
        
        // Solve DPR1 problem
        Freccia::DPR1::DPR1EigenSolver eigh(ew, z, rho);

        // Accumulate result
        ew = eigh.eigenvalues();
        ev *= eigh.eigenvectors();

    }
}

// Eigensystem update problem
void Freccia::DPRK::DPRKEigenSolver::DPRKeigh(const Eigen::Ref<const Eigen::VectorXd> & D, const Eigen::Ref<const Eigen::MatrixXd> & W, const double rho, const Eigen::Ref<const Eigen::MatrixXd> & Q0) {
    // Get dimensions of problem
    assert(D.size() == W.rows());
    unsigned int N = W.rows();
    unsigned int f = W.cols();

    // Initialize result
    ew = D;
    ev = Q0;

    // Compute DPRK update
    for(int i = 0; i < f; ++i){
        // Compute i-th DPR1 update
        Eigen::VectorXd z = ev.transpose()*W.col(i);
        
        // Solve DPR1 problem
        Freccia::DPR1::DPR1EigenSolver eigh(ew, z, rho);

        // Accumulate result
        ew = eigh.eigenvalues();
        ev *= eigh.eigenvectors();
    }
}