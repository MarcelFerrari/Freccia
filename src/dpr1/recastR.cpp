// Eigen
#include <Eigen/Dense>

// Matrix datatypes
#include "Freccia/matrix/matrix.hpp"

// DPR1Solver header
#include "Freccia/dpr1/DPR1EigenSolver.hpp"

void Freccia::DPR1::DPR1EigenSolver::recastR(){
    if(!isRecasted){
        // Recast R to __float128 datatype
        R_ld.D = R.D.cast<__float128>();
        R_ld.z = R.z.cast<__float128>();
        R_ld.rho = static_cast<__float128>(R.rho);
        Rzsqr_ld = R_ld.z * R_ld.z; // Recompute with __float128 precision
        isRecasted = true;
    }
}