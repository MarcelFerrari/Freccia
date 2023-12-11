#ifndef FRECCIA_DPR1EIGENSOLVER_H
#define FRECCIA_DPR1EIGENSOLVER_H

// Standard lib
#include <tuple>

// Eigen
#include <Eigen/Dense>

// Freccia
#include "Freccia/matrix/matrix.hpp"
#include "Freccia/utils/utils.hpp"
#include "Freccia/options/options.hpp"
#include "Freccia/deflation/type1.hpp"
#include "Freccia/deflation/type2.hpp"

namespace Freccia::Arrowhead {

/**
 * Class ArrowheadEigenSolver
 * 
 * Computes all eigenvalues and eigenvectors of a symmetric Arrowhead matrix.
 *
 * This class computes all eigenvalues and eigenvectors of a symmetric Arrowhead matrix in O(n^2)
 * as described in https://doi.org/10.1016/j.laa.2013.10.007. The algorithm applies a preprocessing
 * step known as deflation to reduce the problem to a smaller irreducible one. The deflation step
 * ensures that the elements of D are ordered decreasingly and unique, and that the elements of z
 * are non-zero.
 *
 * @param [Eigen::ArrayXd] D diagonal of the Arrowhead matrix.
 * @param [Eigen::ArrayXd] z shaft of the Arrowhead matrix.
 * @param [double] alpha lower right-corner scalar tip of the Arrowhead shaft.
 * @param [Freccia::Options::ArrowheadEigenSolverOptions] OPTIONS_IN options for the eigensolver.
 */

class ArrowheadEigenSolver {
    // Type aliases
    using Options = Freccia::Options::ArrowheadEigenSolverOptions;
    
    public:
        // Class constructor.
        // If instances of Eigen::VectorXd are passed, they are converted to Eigen::ArrayXd implicitly by Eigen.
        ArrowheadEigenSolver(
            const Eigen::ArrayXd& D, 
            const Eigen::ArrayXd& z,
            double alpha,
            const Options& OPTIONS_IN = Options()
            ) 
        :
        opt(OPTIONS_IN)
        {
            // Search environment variables if necessary
            opt.loadEnv();
            
            // Call eigensolver
            eigh(D, z, alpha);
        }

        // Getter functions
        // Return eigenvalues
        const Eigen::VectorXd& eigenvalues() const { return ew; }
        
        // Return eigenvectors
        const Eigen::MatrixXd& eigenvectors() const { return ev; }

    private:
        // Solver options
        Options opt;
        
        // Problem size
        unsigned int NS; // Size of the source Arrowhead matrix
        unsigned int NR; // Size of the reduced Arrowhead matrix

        // Arrowhead matrices
        ArrowheadMatrix<double> S; // Source Arrowhead matrix
        ArrowheadMatrix<double> R; // Reduced Arrowhead matrix

        // Store R.w^2 to avoid recomputing it at each iteration
        Eigen::ArrayXd Rwsqr; // Square of R.w

        // Quad precision objects
        // The computation of these objects is only necessary for ill-conditioned problems.
        // Quad precision datatypes are slow and do not allow for vectorization, 
        // so recomputing them at each iteration is avoided.
        bool isRecasted = false; // Flag to indicate if R has been recasted to __float128
        ArrowheadMatrix<__float128> R_ld; // Reduced Arrowhead matrix in __float128 datatype
        Eigen::Array<__float128, Eigen::Dynamic, 1> Rwsqr_ld; // __float128 version of Rwsqr

        // Permutation matrix used to sort D and z
        Permutation sort;

        // Type 1 deflation object
        Freccia::Deflation::Type1 type1_deflation;
        
        // Type 2 deflation object
        Freccia::Deflation::Type2 type2_deflation;
        
        // Reduced solution -> this will be optimized and removed in the future
        Eigen::VectorXd DR;
        Eigen::MatrixXd QR;

        // Solution
        Eigen::VectorXd ew; // Eigenwerte (eigenvalues)
        Eigen::MatrixXd ev; // Eigenvektoren (eigenvectors)

        // Internal functions (see source files for documentation)
        void eigh(const Eigen::ArrayXd& D, const Eigen::ArrayXd& z, double alpha);
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

} // namespace Arrowhead
#endif