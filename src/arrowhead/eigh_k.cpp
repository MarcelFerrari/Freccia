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

#include <iomanip> //

// Eigen
#include <Eigen/Dense>

// Matrix datatypes
#include "Freccia/matrix/matrix.hpp"
#include "Freccia/utils/numeric.hpp"

// DPR1Solver header
#include "Freccia/arrowhead/ArrowheadEigenSolver.hpp"

// Shortcut for Freccia::Numeric::sign()
using namespace Freccia::Numeric;

void Freccia::Arrowhead::ArrowheadEigenSolver::eigh_k(unsigned int k){
    // Variables for the sigma (shift), side (left = 0 or right = 1), and index to use in the computation
    bool side;
    double sigma;
    unsigned int i;

    // Depending on the value of k, calculate sigma, side, and i
    if (k == 0){ // First eigenpair
        sigma = R.D(0);  // The shift is the first element of D
        side = true;   // For the first eigenvalue, choose the right side
        i = 0;         // Shift index is 0
    } else if (k == NR - 1) { // Last eigenpair
        sigma = R.D(NR - 2);
        i = NR - 2;
        side = false;
    } else { // Interior eigenvalues, k = 1, ..., NR - 2
        double dk = R.D(k);
        
        // Subtract dk from the diagonal of R
        Eigen::ArrayXd Dtmp = R.D - dk; // Avoid potential aliasing
        double btmp = R.b - dk;

        // Compute one iteration of bisection to find the shift side and index
        double middle = Dtmp(k-1)/2.;    // Compute the shift tau

        // Execute one iteration of bisection to find the shift side and index
        double F = btmp - middle - (Rwsqr/(Dtmp - middle)).sum();

        // Depending on the value of F, calculate sigma, side, and i
        if(F < 0.){
            sigma = R.D(k);
            side = true;
            i = k;
        } else {
            sigma = R.D(k-1);
            side = false;
            i = k - 1;
        }
    }

    // Compute shifted inverse of R with shift sigma = R.D(i)
    ArrowheadMatrix<double> Rinv = shiftInvert(i);

    // Compute nu
    double nu = solveSecularEQ(Rinv, side);

    
    // Compute Knu and nu1
    auto [Knu, nu1] = computeKnu(Rinv, nu);

    // Check if we have to recompute the eigenvalue
    unsigned int niter = 0;
    while(Knu > opt.KNU_TOL && niter < opt.RECOMPUTE_MAX_ITER){ // Knu >> 1
        //Shift between lambda and pole
        std::tuple<double, double, double> rax = recomputeEigenvalue(nu, nu1, sigma, side);
        nu = std::get<0>(rax);
        nu1 = std::get<1>(rax);
        sigma = std::get<2>(rax);
        
        // Recompute Knu
        Knu = nu1 / std::abs(nu);
        ++niter;
    }
    
    // Calculate mu and the eigenvalue lambda
    double mu = 1./nu;
    double lambda = mu + sigma;

    // Compute the corresponding eigenvector
    vect(sigma, mu, k);
    
    // Check if lambda needs to be recomputed
    if((std::abs(R.D(i)) + std::abs(mu)) / std::abs(lambda) > opt.K_LAMBDA_TOL){
       if ( (k == 0 && R.D[0] < 0.0) || 
            (k == NR - 1 && R.D[NR - 2] > 0.0) || 
            (i < NR - 2 && !side && (sign(R.D[i]) + sign(R.D[i + 1]) == 0)) || 
            (i > 0 && side && (sign(R.D[i]) + sign(R.D[i - 1]) == 0)))
        {
            lambda = recomputeEigenvalue(k);
        }
    }
    
    
    DR(k) = lambda;
}