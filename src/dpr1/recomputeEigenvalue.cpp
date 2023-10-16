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

// DPR1Solver header
#include "Freccia/dpr1/DPR1EigenSolver.hpp"

std::pair<double, double> Freccia::DPR1::DPR1EigenSolver::computeKnu(const ArrowheadMatrix<double> &Rinv, const double nu){
    // Compute Knu
    double nu1 = 0.0;
    
    // Compute nu1 = ||(D + zz^T - sigma I)^-1||_1,2
    for (unsigned int j = 0; j < NR - 1; j++) {
        if(j == Rinv.i){continue;} // Skip the index i
        nu1 = std::max(nu1, std::abs(Rinv.D(j)) + std::abs(Rinv.w(j)));
    }

    nu1 = std::max(nu1, Rinv.w.abs().sum() + std::abs(Rinv.b));
    double Knu = nu1 / std::abs(nu);

    return std::make_pair(Knu, nu1);
}

std::tuple<double, double, double> Freccia::DPR1::DPR1EigenSolver::recomputeEigenvalue(double nu, double nu1, double sigma_zero, bool side){
    // Compute new shift with magic
    nu = (side == true) ? std::abs(nu) : -std::abs(nu);
    nu1 = (nu > 0.) ? -nu1 : nu1;
    double sigma = (nu + nu1)/(2.0 * nu * nu1) + sigma_zero;

    // Compute new Rinv
    DPR1Matrix<double> Rinv = shiftInvert(sigma); // This time Rinv is a full size DPR1 Matrix
    
    // Compute bounds for the eigenvalue
    double left, right;
    Eigen::ArrayXd zabs = Rinv.z.abs();
    Eigen::ArrayXd zsqr = Rinv.z * Rinv.z;
    
    // Sort the D array decreasingly
    Permutation p;
    p.argsort(Rinv.D);
    const std::vector<unsigned int> & idx = p.perm();

    if(Rinv.rho > 0.0){ // Standard case
        if (side == true){ // Right side
            left = Rinv.D(idx[0]);
            right = Rinv.D(idx[0]) + Rinv.rho * zsqr.sum();
        } else { // Left side
            left = Rinv.D(idx[NR-1]);
            right = Rinv.D(idx[NR-2]);
        }
    } else { // rho < 0 case
        if(side == true){ // Right side
            left = Rinv.D(idx[1]);
            right = Rinv.D(idx[0]);
        } else { // Left side
            left = Rinv.D(idx[NR-1]) + Rinv.rho * zsqr.sum();
            right = Rinv.D(idx[NR-1]);
        }
    }
    
    // Define the secular function to be solved
    auto F = [&](double x){return 1. + Rinv.rho * (zsqr/(Rinv.D - x)).sum();};

    double middle = (left + right) / 2.;
    unsigned niter = 0;
    double eps = std::numeric_limits<double>::epsilon();
    double sgn = (Rinv.rho > 0.0) ? 1.0 : -1.0;
    
    // Recompute nu1
    nu1 =  Rinv.D.abs().maxCoeff() + std::abs(Rinv.rho)*zsqr.sum();

    // Recompute nu
    while(((right - left) > 2.*eps*std::max(std::abs(left), std::abs(right))) && niter < opt.BISECT_MAX_ITER){
        
        if(sgn*F(middle) < 0.){
            left = middle;
        } else {
            right = middle;
        }
        
        middle = (left + right) / 2.;
        ++niter;
    }

    return std::make_tuple(right, nu1, sigma);    
};

double Freccia::DPR1::DPR1EigenSolver::recomputeEigenvalue(unsigned int k){
    // Compute Lambda via rayleigh quotient using the computed eigenvector q
        DPR1Matrix<double> Rinv = shiftInvert(0.0);
        
        if(std::isinf(Rinv.rho)){ // Check if matrix is singular
            return 0.0; 
        } else { // Compute lambda as Rayleigh quotient
            // qT(D + rho * vvT)q = qTDq + rho * (qTz) (zTq)
            auto view = ev(type2_deflation.nnzero(), type2_deflation.nnzero());
            double qTz = view.col(k).dot(Rinv.z.matrix());
            double qTDq = view.col(k).dot(Rinv.D.matrix().cwiseProduct(view.col(k)));
        
            return 1./(qTDq + Rinv.rho*qTz*qTz);
        }
}