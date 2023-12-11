#ifndef FRECCIA_DPR1EIGENSOLVER_H
#define FRECCIA_DPR1EIGENSOLVER_H

// Standard library includes
#include <tuple>

// Eigen library for linear algebra
#include <Eigen/Dense>

// Freccia includes
#include "Freccia/matrix/matrix.hpp"
#include "Freccia/utils/utils.hpp"
#include "Freccia/options/options.hpp"
#include "Freccia/deflation/type1.hpp"
#include "Freccia/deflation/type2.hpp"

namespace Freccia::DPR1 {

/**
 * Class DPR1EigenSolver
 * 
 * Computes all eigenvalues and eigenvectors of a symmetric DPR1 (Diagonal Plus Rank-One) matrix.
 *
 * This class computes all eigenvalues and eigenvectors of a symmetric DPR1 matrix in O(n^2)
 * as described in https://doi.org/10.1016/j.laa.2015.09.025. The algorithm applies a preprocessing
 * step known as deflation to reduce the problem to a smaller irreducible one. The deflation step
 * ensures that the elements of D are ordered decreasingly and unique, and that the elements of z
 * are non-zero.
 *
 * @param [Eigen::ArrayXd] D diagonal elements of the DPR1 matrix.
 * @param [Eigen::ArrayXd] z elements of rank one update vector of the DPR1 matrix.
 * @param [double] rho scalar multiplier for the rank-one update.
 * @param [Freccia::Options::DPR1EigenSolverOptions] OPTIONS_IN options for the eigensolver.
 */
class DPR1EigenSolver {
    // Type aliases for readability
    using Options = Freccia::Options::DPR1EigenSolverOptions;
    
public:
    // Constructor
    // Implicit conversion from Eigen::VectorXd to Eigen::ArrayXd is supported by Eigen.
    DPR1EigenSolver(
        const Eigen::ArrayXd& D, 
        const Eigen::ArrayXd& z,
        double rho = 1.0,
        const Options& OPTIONS_IN = Options()
                    ) 

    : opt(OPTIONS_IN)
    {
        // Load environment variables if necessary
        opt.loadEnv(); 
        
        // Call the eigensolver driver function
        eigh(D, z, rho);
    }

    // Getter functions
    // Returns eigenvalues of the DPR1 matrix
    const Eigen::VectorXd& eigenvalues() const { return ew; }
    
    // Returns eigenvectors of the DPR1 matrix
    const Eigen::MatrixXd& eigenvectors() const { return ev; }

private:
    // Solver options
    Options opt;
    
    // Problem size
    unsigned int NS; // Size of the source DPR1 matrix
    unsigned int NR; // Size of the reduced DPR1 matrix

    // DPR1 matrices
    DPR1Matrix<double> S; // Source DPR1 matrix
    DPR1Matrix<double> R; // Reduced DPR1 matrix
    Eigen::ArrayXd Rzsqr; // Square of R.z for efficient computation

    // Quad precision objects for ill-conditioned problems
    bool isRecasted = false; // Flag for recasting to __float128
    DPR1Matrix<__float128> R_ld; // Reduced DPR1 matrix in __float128 datatype
    Eigen::Array<__float128, Eigen::Dynamic, 1> Rzsqr_ld; // __float128 version of Rzsqr

    // Permutation matrix for sorting D and z
    Permutation sort;

    // Deflation objects
    Freccia::Deflation::Type1 type1_deflation;
    Freccia::Deflation::Type2 type2_deflation;
    
    // Solution storage
    Eigen::VectorXd ew; // Eigenwerte (eigenvalues)
    Eigen::MatrixXd ev; // Eigenvektoren (eigenvectors)

    // Internal functions (see source files for documentation)
    void eigh(const Eigen::ArrayXd& D, const Eigen::ArrayXd& z, double rho);
    void eigh_k(unsigned int k);
    void recastR();
    ArrowheadMatrix<double> shiftInvert(unsigned int i);
    DPR1Matrix<double> shiftInvert(double sigma);
    double solveSecularEQ(const ArrowheadMatrix<double>& Rinv, bool side);
    void vect(double sigma, double mu, unsigned int k);
    std::pair<double, double> computeKnu(const ArrowheadMatrix<double> &Rinv, const double nu);
    std::tuple<double, double, double> recomputeEigenvalue(double nu, double nu1, double sigma_zero, bool side);
    double recomputeEigenvalue(unsigned int k);
};

} // namespace Freccia::DPR1
#endif
