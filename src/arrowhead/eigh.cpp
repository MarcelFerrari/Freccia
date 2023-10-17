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
#include "Freccia/arrowhead/ArrowheadEigenSolver.hpp"

static void printMat(ArrowheadMatrix<double>& M){
    std::cout << "D = " << M.D.transpose() << std::endl;
    std::cout << "w = " << M.w.transpose() << std::endl;
    std::cout << "b = " << M.b << std::endl;
}

// Main driver function
void Freccia::Arrowhead::ArrowheadEigenSolver::eigh(const Eigen::ArrayXd& D, const Eigen::ArrayXd& z, double alpha) {        
    
    // Quick return if z vector is zero
    if(z.matrix().isZero()){
        // Init eigenvalues and eigenvectors
        ew = Eigen::VectorXd::Zero(D.size()+1);
        ev = Eigen::MatrixXd::Identity(D.size(), D.size());

        // Copy D and alpha into ew
        ew.head(D.size()) = D;
        ew(D.size()) = alpha;
        
        return;
    }
    
    // Sort input
    sort.argsort(D);

    // Initialize the source matrix
    S.D = D(sort.perm());
    S.w = z(sort.perm());
    S.b = alpha;

    // Source matrix size
    NS = D.size() + 1;

    // Initialize the eigenvalues and eigenvectors
    ew = Eigen::VectorXd::Zero(NS);
    ev = Eigen::MatrixXd::Zero(NS, NS);

    // Compute and apply type 1 deflation to S
    type1_deflation.deflate(S.D, S.w, opt);

    // Compute type 2 deflation of S
    type2_deflation.deflate(S.D, S.w, ew, ev);
    
    // Build reduced arrowhead matrix R
    // Recasting R to __float128 will be done lazily if necessary
    const std::vector<unsigned int>& nnzero = type2_deflation.getPartition().nnzero(); // Final non-zero elements indices
    NR = nnzero.size() + 1;

    // Set up reduced solution
    DR = Eigen::VectorXd::Zero(NR);
    QR = Eigen::MatrixXd::Zero(NR, NR);
    
    // Quick return if there is one non-zero element
    if(NR == 2){ // 2 because NR is the full size of the arrowhead matrix including b
        // Construct 2x2 dense matrix representation of the arrowhead matrix
        Eigen::Matrix<double, 2, 2> M;
    
        M << S.D(nnzero[0]), S.w(nnzero[0]),
             S.w(nnzero[0]), S.b;

        // Use Eigen's selfadjointeigensolver that relies on closed form solution
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M);

        // Copy eigenvalues and eigenvectors
        ew(nnzero) = solver.eigenvalues();
        ev(nnzero, nnzero) = solver.eigenvectors();

    } else { // Standard computation
        // Create reduced problem matrix for convenience        
        R.D = S.D(nnzero);
        R.w = S.w(nnzero);
        R.b = S.b;
        Rwsqr = R.w * R.w; // Useful for later

        // Compute the eigenpairs of the reduced problem
        #pragma omp parallel for schedule(dynamic)
        for(unsigned int i = 0; i < NR; ++i){
            eigh_k(i);
        }
    }

    // Copy solution back to ew and ev
    ew(nnzero) = DR;
    ew(NS - 1) = DR(NR - 1); // Last element
    
    // Not so trivial
    ev(nnzero, nnzero) = QR.topLeftCorner(NR - 2, NR - 2); // Top left corner
    ev.col(NS-1)(nnzero) = QR.col(NR - 1).head(NR - 1); // Last column
    ev.row(NS-1)(nnzero) = QR.row(NR - 1).head(NR - 1); // Last row
    ev(NS - 1, NS - 1) = QR(NR - 1, NR - 1); // Last element
    
    // Apply type 1 deflation to the eigenvectors
    // The HH matrix is templated so we can pass Eigen View object directly
    // This is a hack and should be fixed in the future
    const HHDeflationMatrix& type1_HH = type1_deflation.getHH();
    const std::vector<unsigned int>& type1_nnzero = type1_deflation.getPartition().nnzero();
    
    // Apply to nnzero block
    type1_HH.applyToTheLeft(ev(type1_nnzero, type1_nnzero));
    
    // ONLY apply to last column of ev
    Eigen::VectorXd q = ev.col(NS - 1)(type1_nnzero);
    type1_HH.applyToVec(q);
    ev.col(NS - 1)(type1_nnzero) = q;

    // Finally, unsort the eigenvalues and eigenvectors
    const std::vector<unsigned int> & inv = sort.invperm();
    
    // First deal with eigenvalues
    ew.head(NS - 1) = ew.head(NS - 1)(inv).eval();
    
    // Then deal with eigenvectors
    ev.topLeftCorner(NS - 1, NS - 1) = ev.topLeftCorner(NS - 1, NS - 1)(inv, inv).eval(); // Top left corner
    ev.col(NS-1).head(NS-1) = ev.col(NS-1).head(NS-1)(inv).eval(); // Last column
    ev.row(NS-1).head(NS-1) = ev.row(NS-1).head(NS-1)(inv).eval(); // Last row
    
}