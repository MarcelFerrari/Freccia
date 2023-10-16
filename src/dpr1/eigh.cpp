// Standard library
#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <string>
#include <tuple>

// Eigen
#include <Eigen/Dense>

// Matrix datatypes
#include "Freccia/matrix/matrix.hpp"

// DPR1Solver header
#include "Freccia/dpr1/DPR1EigenSolver.hpp"

// Main driver function
void Freccia::DPR1::DPR1EigenSolver::eigh(const Eigen::ArrayXd& D, const Eigen::ArrayXd& z, double rho) {        
    
    // Quick return if z vector is zero
    if(z.matrix().isZero()){
        ew = D;
        ev = Eigen::MatrixXd::Identity(D.size(), D.size());
        return;
    }

    // Check for rho < 0
    bool flip = rho < 0.0;

    // Store the input parameters and flip D if needed
    if(flip){
        S.D = -D;
        S.z = z;
        S.rho = -rho;
    } else {
        S.D = D;
        S.z = z;
        S.rho = rho;
    }

    // Source matrix size
    NS = D.size();

    // Initialize the eigenvalues and eigenvectors
    ew = Eigen::VectorXd::Zero(NS);
    ev = Eigen::MatrixXd::Zero(NS, NS);

    // Sort the source matrix
    sort.argsort(S.D);

    // Apply permutation to S
    S.D = S.D(sort.perm()).eval();
    S.z = S.z(sort.perm()).eval();
    
    // Compute and apply type 1 deflation to S
    type1_deflation.deflate(S.D, S.z, opt);

    // Compute type 2 deflation of S
    type2Deflation();
    
    // Build reduced DPR1 matrix R
    // Recasting R to __float128 will be done lazily if necessary
    const std::vector<unsigned int>& nnzero = type2_deflation.nnzero(); // Final non-zero elements indices
    NR = nnzero.size();
    
    // Quick return if there is one non-zero element
    if(NR == 1){
        if(std::abs(S.D(nnzero[0])) < 1e-15 && (S.rho < 1e-15 || std::abs(S.z(nnzero[0])) < 1e-15)){
            ew(nnzero[0]) = 0.;
        } else {
            ew(nnzero[0]) = S.D(nnzero[0]) + S.rho * S.z(nnzero[0]) * S.z(nnzero[0]);
        }
        
        ev(nnzero[0], nnzero[0]) = 1.;
    } else { // Standard computation
        // Create reduced problem matrix for convenience        
        R.D = S.D(nnzero);
        R.z = S.z(nnzero);
        R.rho = S.rho;
        Rzsqr = R.z * R.z; // Useful for later

        // Compute the eigenpairs of the reduced problem
        #pragma omp parallel for schedule(dynamic)
        for(unsigned int i = 0; i < NR; ++i){
            eigh_k(i);
        }
    }
    
    // Apply type 1 deflation to the eigenvectors
    // Easier to make a copy of the eigenvectors and then apply the deflation
    // Only reconstruct Q if needed!
    if(opt.RECONSTRUCT_Q){
        // The HH matrix takes an Eigen::Ref object so we can pass an IndexedView directly
        { // Do not pollute namespace
            const std::vector<unsigned int>& type1_nnzero = type1_deflation.getPartition().nnzero();
            const HHDeflationMatrix& type1_HH = type1_deflation.getHH();
            type1_HH.applyToTheLeft(ev(type1_nnzero, type1_nnzero));
        }
    }

    // Finally, unsort the eigenvalues and eigenvectors
    ew = (flip) ? -ew(sort.invperm()).eval() : ew(sort.invperm()).eval();
    ev = ev(sort.invperm(), sort.invperm()).eval();
}