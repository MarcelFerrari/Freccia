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
#include "Freccia/utils/root_finding.hpp"

double Freccia::DPR1::DPR1EigenSolver::solveSecularEQ(const ArrowheadMatrix<double>& Rinv, bool side){
    // Initialize the left and right bounds of the solution
    double left, right;
    
    // Compute the absolute value of each element in vector w
    Eigen::ArrayXd wabs = Rinv.w.abs();
   
    // Depending on the side, compute the initial bounds for the eigenvalue
    if(side == false){ // Left side
        left = std::min((Rinv.D-wabs).minCoeff(), Rinv.b - wabs.sum());
        right = Rinv.D.minCoeff();
    } else { // Right side
        right = std::max((Rinv.D+wabs).maxCoeff(), Rinv.b + wabs.sum());
        left = Rinv.D.maxCoeff();
    }
    
    // Compute the square of each element in vector w
    Eigen::ArrayXd wsqr = Rinv.w * Rinv.w;

    // Define the secular function to be solved
    auto F = [&](double x){return Rinv.b - x - (wsqr / (Rinv.D - x)).sum();};

    // Try with fast solver first
    {
        double nu = Freccia::RootFinding::toms748(left, right, F, opt);

        if(!std::isnan(nu)){ // Check that the fast solver did not fail
            return nu;
        }
    }

    // Fast solver failed, resort to bisection
    double middle = (left + right) / 2.;
    unsigned niter = 0;
    double eps = std::numeric_limits<double>::epsilon();
    
    while((right - left) > 2.*eps*std::max(std::abs(left), std::abs(right)) && niter < opt.BISECT_MAX_ITER){
        if(F(middle) > 0){
            left = middle;
        } else {
            right = middle;
        }
        middle = (left + right) / 2.;
        ++niter;
    }

    return right;
}