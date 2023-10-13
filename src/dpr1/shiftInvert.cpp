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

ArrowheadMatrix<double> Freccia::DPR1::DPR1EigenSolver::shiftInvert(unsigned int i){
    // Compute the shift sigma
    double sigma = R.D(i);
    
    // Compute shifted matrix Dshift
    Eigen::ArrayXd Dshift = R.D - sigma;

    // Compute the inverse of the ith element of vector z
    double w_xi = 1./R.z(i);
    
    // Check conditioning of problem
    // Compute Kz
    double Kz = (R.z.segment(0, i).abs().sum() + R.z.segment(i+1, NR - i - 1).abs().sum())*w_xi;
    
    // Compute Kb
    // The inverse is partitioned into two segments
    double P = 0.;
    double Q = 0.;
    
    ArrowheadMatrix<double> Rinv(NR-1);
    Rinv.i = i;

    // Compute invD1 and zinvD1 if the first partition is non-empty
    if(i > 0){
        Rinv.D.segment(0, i) = 1./Dshift.segment(0,i); // First segment of the inverse of the shifted matrix
        P = (Rzsqr.segment(0, i) * Rinv.D.segment(0, i)).sum(); // Auxiliary variable for the computation of b
        Rinv.w.segment(0, i) = - w_xi * (Rinv.D.segment(0, i) * R.z.segment(0, i)); // First segment of the computed vector w
    }

    // Compute invD2, w2 and zinvD2 if i is less than N-1 if the second partition is non-empty
    if(i < NR - 1){
        Rinv.D.segment(i, NR - i - 1) = 1./Dshift.segment(i+1, NR - i - 1); // Second segment of the inverse of the shifted matrix
        Q = (Rzsqr.segment(i+1, NR - i - 1) * Rinv.D.segment(i, NR - i - 1)).sum(); // Auxiliary variable for the computation of b
        Rinv.w.segment(i, NR - i - 1) = - w_xi * (Rinv.D.segment(i, NR - i - 1) * R.z.segment(i+1, NR - i - 1)); // Second segment of the computed vector w
    }
    
    R.rho > 0.0 ? P += 1./R.rho : Q += 1./R.rho; // Deal with last element

    // Compute Kb
    double Kb = (P - Q)/std::abs(P + Q);

    // Standard precision is enough
    if(Kz < opt.KZ_TOL*NR || Kb < opt.KB_TOL){
        // Compute b using the previously calculated values
        Rinv.b = (P+Q)*w_xi*w_xi;
    } else { // Need to recompute b using quad precision
        // Recast R to __float128 if needed
        recastR();

        // Recast input to __float128 datatype and recompute b
        Eigen::Array<__float128, Eigen::Dynamic, 1> Dshift_ld = R_ld.D - static_cast<__float128>(sigma);

        __float128 P_ld = 0.0L;
        __float128 Q_ld = 0.0L;
        
        // Compute invD1_ld and P_ld if the first partition is non-empty
        if(i > 0){
            P_ld = (Rzsqr_ld.segment(0, i) / Dshift_ld.segment(0, i)).sum(); // Auxiliary variable for the computation of b
        }

        // Compute invD2_ld, w2 and Q_ld if i is less than N-1 if the second partition is non-empty
        if(i < NR-1){
            Q_ld = (Rzsqr_ld.segment(i+1, NR - i - 1) / Dshift_ld.segment(i+1, NR - i - 1)).sum(); // Auxiliary variable for the computation of b
        }

        R.rho > 0.0 ? P_ld += 1.0L/R_ld.rho : Q_ld += 1.0L/R_ld.rho; // Deal with last element
        Rinv.b = static_cast<double>((P_ld + Q_ld)/(Rzsqr_ld(i)));
    }

    return std::move(Rinv);
}

DPR1Matrix<double> Freccia::DPR1::DPR1EigenSolver::shiftInvert(double sigma){
    DPR1Matrix<double> Rinv;
    // Try double precision
    Rinv.D = 1./(R.D - sigma);
    Rinv.z = Rinv.D * R.z;

    double P = 0.0, Q = 0.0;
    Eigen::ArrayXd zsqrDinv = Rzsqr * Rinv.D;
    for(unsigned int i = 0; i < NR; ++i){
        if(zsqrDinv(i) > 0.0){
            P += zsqrDinv(i);
        } else {
            Q += zsqrDinv(i);
        }
    }

    // Deal with rho
    R.rho > 0.0 ? P += 1/R.rho : Q += 1/R.rho;

    // Check conditioning of rho
    double Krho = (P-Q)/std::abs(P+Q);
    if(Krho < opt.RHO_TOL){ // Double precision is accurate enough
        Rinv.rho = -1. / (P+Q);
    } else { // Quad precision is required
        // Recast R to __float128 if needed
        recastR();

        // We only need to store z_ld
        Eigen::Array<__float128, Eigen::Dynamic, 1> z_ld = Rzsqr_ld / (R_ld.D - static_cast<__float128>(sigma));

        // Explicit loop here as it is faster than using unaryExpr
        __float128 Pd = 0.0L, Qd = 0.0L;

        for(unsigned int i = 0; i < NR; ++i){
            if(z_ld(i) > 0.0L){
                Pd += z_ld(i);
            } else {
                Qd += z_ld(i);
            }
        }

        R.rho > 0.0 ? Pd += 1.0L/R_ld.rho : Qd += 1.0L/R_ld.rho;
        Rinv.rho = static_cast<double>(-1.0L/(Pd+Qd));
    }

    return std::move(Rinv);
}