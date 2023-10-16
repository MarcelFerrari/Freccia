// Eigen
#include <Eigen/Dense>

// Matrix datatypes
#include "Freccia/matrix/matrix.hpp"

// DPR1Solver header
#include "Freccia/dpr1/DPR1EigenSolver.hpp"


void Freccia::DPR1::DPR1EigenSolver::vect(double sigma, double mu, unsigned int k){
    
    // Get view to eigenvectors
    const std::vector<unsigned int>& perm = type2_deflation.getPartition().nnzero();
    auto view = ev(perm, perm);
    
    // Compute eigenvector
    view.col(k) = (R.z / ((R.D - sigma) - mu)).matrix();
    
    // Normalize the resulting vector to obtain v = x / ||x||_2
    view.col(k).normalize();
}