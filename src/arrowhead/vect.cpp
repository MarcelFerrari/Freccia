// Eigen
#include <Eigen/Dense>

// Matrix datatypes
#include "Freccia/matrix/matrix.hpp"

// DPR1Solver header
#include "Freccia/arrowhead/ArrowheadEigenSolver.hpp"

void Freccia::Arrowhead::ArrowheadEigenSolver::vect(double sigma, double mu, unsigned int k){
    // Compute eigenvector
    QR.col(k).head(NR-1) = (R.w / ((R.D - sigma) - mu)).matrix();
    QR(NR-1, k) = - 1.;
    
    // Normalize the resulting vector to obtain v = x / ||x||_2
    QR.col(k).normalize();
}