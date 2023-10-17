// Eigen
#include <Eigen/Dense>

// Matrix datatypes
#include "Freccia/matrix/matrix.hpp"

// DPR1Solver header
#include "Freccia/arrowhead/ArrowheadEigenSolver.hpp"

void Freccia::Arrowhead::ArrowheadEigenSolver::recastR(){
    if(!isRecasted){
        // Recast R to __float128 datatype
        R_ld.D = R.D.cast<__float128>();
        R_ld.w = R.w.cast<__float128>();
        R_ld.b = static_cast<__float128>(R.b);
        Rwsqr_ld = R_ld.w * R_ld.w; // Recompute with __float128 precision
        isRecasted = true;
    }
}