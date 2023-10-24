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

class ArrowheadEigenSolver {
    public:
        // Constructor for ArrayXd
        ArrowheadEigenSolver(const Eigen::ArrayXd& D, 
                        const Eigen::ArrayXd& z,
                        double alpha,
                        const Freccia::Options::ArrowheadEigenSolverOptions& OPTIONS_IN
                                = Freccia::Options::ArrowheadEigenSolverOptions()
                        ) 
        :
        opt(OPTIONS_IN)
        {
            // Search environment variables if necessary
            opt.loadEnv();
            // Call eigensolver after setting up tolerances
            eigh(D, z, alpha);
        }

        // Return eigenvalues
        const Eigen::VectorXd& eigenvalues() const { return ew; }
        
        // Return eigenvectors
        const Eigen::MatrixXd& eigenvectors() const { return ev; }

    private:
        // Solver options
        Freccia::Options::ArrowheadEigenSolverOptions opt;
        
        // Input parameters
        int NS;
        int NR;

        ArrowheadMatrix<double> S; // Source DPR1 matrix
        
        ArrowheadMatrix<double> R; // Reduced DPR1 matrix
        Eigen::ArrayXd Rwsqr; // Square of R.w

        bool isRecasted = false; // Flag to indicate if R has been recasted to __float128
        ArrowheadMatrix<__float128> R_ld; // Reduced DPR1 matrix in __float128 datatype
        Eigen::Array<__float128, Eigen::Dynamic, 1> Rwsqr_ld; // __float128 version of Rwsqr

        // Sorting
        Permutation sort;

        // Type 1 deflation
        Freccia::Deflation::Type1 type1_deflation;
        
        // Type 2 deflation
        Freccia::Deflation::Type2 type2_deflation;
        
        // Reduced solution
        Eigen::VectorXd DR;
        Eigen::MatrixXd QR;

        // Output parameters
        Eigen::VectorXd ew; // Eigenwerte
        Eigen::MatrixXd ev; // Eigenvektoren

        // Functions
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

} // namespace DPR1
#endif