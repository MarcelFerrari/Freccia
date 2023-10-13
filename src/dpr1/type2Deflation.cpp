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

void Freccia::DPR1::DPR1EigenSolver::type2Deflation(){
    // Partition zero and non zero values in S.z
    type2_deflation.partition(S.z);
    
    // Compute the deflated solution. For each index 'k' corresponding to a zero 
    // element in 'z', set the (k,k)-th element of matrix Q to 1 and the k-th 
    // element of vector L to the corresponding element in D_sorted. This effectively 
    // sets the diagonal elements of Q and L for zero elements, and leaves non-zero 
    // elements to be handled separately.
    for(unsigned int k: type2_deflation.zero()){
        ev(k, k) = 1.;
        ew(k) = S.D(k);
    }
}