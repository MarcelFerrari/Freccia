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
#include "Freccia/arrowhead/ArrowheadEigenSolver.hpp"

ArrowheadMatrix<double> Freccia::Arrowhead::ArrowheadEigenSolver::shiftInvert(unsigned int i){
    // Compute the shift sigma
    double sigma = R.D(i);
    
    // Compute shifted matrix Dshift
    Eigen::ArrayXd Dshift = R.D - sigma;
    double bshift = R.b - sigma;

    // Compute the inverse of the ith element of vector z
    double w_xi = 1./R.w(i);
    
    // Check conditioning of problem
    // Compute Kz
    double Kz = (R.w.segment(0, i).abs().sum() + R.w.segment(i+1, NR - i - 2).abs().sum())*w_xi;
    
    // Compute Kb
    // The inverse is partitioned into two segments
    double P = 0.;
    double Q = 0.;
    
    // NR - 1 as the ArrowheadMatrix constructor takes the size of D and z
    ArrowheadMatrix<double> Rinv(NR-1);
    Rinv.i = i;

    // Compute invD1 and zinvD1 if the first partition is non-empty
    if(i > 0){
        Rinv.D.segment(0, i) = 1./Dshift.segment(0,i); // First segment of the inverse of the shifted matrix
        P = (Rwsqr.segment(0, i) * Rinv.D.segment(0, i)).sum(); // Auxiliary variable for the computation of b
        Rinv.w.segment(0, i) = - w_xi * (Rinv.D.segment(0, i) * R.w.segment(0, i)); // First segment of the computed vector w
    }

    // Compute invD2, w2 and zinvD2 if i is less than N-1 if the second partition is non-empty
    if(i < NR - 1){
        // NR - i - 2 as the ArrowheadMatrix constructor takes the size of D and z and not the whole matrix!
        Rinv.D.segment(i+1, NR - i - 2) = 1./Dshift.segment(i+1, NR - i - 2); // Second segment of the inverse of the shifted matrix
        Q = (Rwsqr.segment(i+1, NR - i - 2) * Rinv.D.segment(i+1, NR - i - 2)).sum(); // Auxiliary variable for the computation of b
        Rinv.w.segment(i+1, NR - i - 2) = - w_xi * (Rinv.D.segment(i+1, NR - i - 2) * R.w.segment(i+1, NR - i - 2)); // Second segment of the computed vector w
    }
    

    Rinv.w(i) = w_xi; // Middle element of w

    // Deal with R.b
    bshift > 0.0 ? Q -= bshift : P -= bshift;
    
    // Compute Kb
    double Kb = (P - Q)/std::abs(P + Q);

    // Standard precision is enough
    if((Kz < opt.KZ_TOL*NR || Kb < opt.KB_TOL)){
        // Compute b using the previously calculated values
        Rinv.b = (P+Q)*w_xi*w_xi;
    } else { // Need to recompute b using quad precision
        // Recast R to __float128 if needed
        recastR();

        // Recast input to __float128 datatype and recompute b
        Eigen::Array<__float128, Eigen::Dynamic, 1> Dshift_ld = R_ld.D - static_cast<__float128>(sigma);

        __float128 bshift_ld = R_ld.b - static_cast<__float128>(sigma);
        __float128 P_ld = 0.0L;
        __float128 Q_ld = 0.0L;
        
        // Compute invD1_ld and P_ld if the first partition is non-empty
        if(i > 0){
            P_ld = (Rwsqr_ld.segment(0, i) / Dshift_ld.segment(0, i)).sum(); // Auxiliary variable for the computation of b
        }

        // Compute invD2_ld, w2 and Q_ld if i is less than N-1 if the second partition is non-empty
        if(i < NR - 1){
            Q_ld = (Rwsqr_ld.segment(i+1, NR - i - 2) / Dshift_ld.segment(i+1, NR - i - 2)).sum(); // Auxiliary variable for the computation of b
        }

        // Deal with R.b
        bshift_ld > 0.0 ? Q_ld -= bshift_ld : P_ld -= bshift_ld;

        Rinv.b = static_cast<double>((P_ld + Q_ld)/(Rwsqr_ld(i)));
    }

    return std::move(Rinv);
}


DPR1Matrix<double> Freccia::Arrowhead::ArrowheadEigenSolver::shiftInvert(double sigma){
    // Inverse is a NR x NR matrix
    DPR1Matrix<double> Rinv(NR);
    
    double bshift = R.b - sigma;

    // Compute inverse
    Rinv.D.head(NR-1) = 1./(R.D - sigma);
    Rinv.z.head(NR-1) = Rinv.D * R.w;
    
    // Compute contribution to P and Q
    Eigen::ArrayXd zsqrDinv = Rwsqr * Rinv.D;
    
    double P = zsqrDinv.unaryExpr([](double x) { return x > 0.0 ? x : 0.0; }).sum();
    double Q = zsqrDinv.unaryExpr([](double x) { return x <= 0.0 ? x : 0.0; }).sum();
    
    // Compute P and Q with explicit for loops
    // for(unsigned int j = 0; j < NR-1; ++j){
    //     if(zsqrDinv(j) > 0.0){
    //         P += zsqrDinv(j);
    //     } else {
    //         Q += zsqrDinv(j);
    //     }
    // }

    Rinv.D(NR-1) = 0.; // i-th element is 0
    Rinv.z(NR-1) = -1.; // i-th element is -1

    // Deal with bshift
    bshift > 0.0 ? Q -= bshift : P -= bshift;

    // Check conditioning of rho
    double Krho = (P-Q)/std::abs(P+Q);

    if(Krho < opt.RHO_TOL){ // Double precision is accurate enough
        Rinv.rho = -1. / (P+Q);
    } else { // Quad precision is required
        // Recast R to __float128 if needed
        recastR();

        // We only need to store w_ld
        Eigen::Array<__float128, Eigen::Dynamic, 1> z_ld = Rwsqr_ld / (R_ld.D - static_cast<__float128>(sigma));
        __float128 bshift_ld = R_ld.b - static_cast<__float128>(sigma);

        // Explicit loop here as it is faster than using unaryExpr as no vectorization is possible
        __float128 P_ld = 0.0L, Q_ld = 0.0L;

        for(unsigned int i = 0; i < NR - 1; ++i){
            if(z_ld(i) > 0.0L){
                P_ld += z_ld(i);
            } else {
                Q_ld += z_ld(i);
            }
        }

        bshift_ld > 0.0 ? Q_ld -= bshift_ld : P_ld -= bshift_ld;
        Rinv.rho = static_cast<double>(-1.0L/(P_ld+Q_ld));
    }

    return std::move(Rinv);
}
