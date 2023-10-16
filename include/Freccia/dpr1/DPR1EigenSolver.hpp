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

namespace Freccia::DPR1 {

class DPR1EigenSolver {
    public:
        // Constructor for ArrayXd
        DPR1EigenSolver(const Eigen::ArrayXd& D, 
                        const Eigen::ArrayXd& z,
                        double rho = 1.0,
                        const Freccia::Options::DPR1EigenSolverOptions& OPTIONS_IN
                                = Freccia::Options::DPR1EigenSolverOptions()
                        ) 
        :
        opt(OPTIONS_IN)
        {
            // Search environment variables if necessary
            opt.loadEnv();
            // Call eigensolver after setting up tolerances
            eigh(D, z, rho);
        }

        // Return eigenvalues
        const Eigen::VectorXd& eigenvalues() const { return ew; }
        
        // Return eigenvectors
        const Eigen::MatrixXd& eigenvectors() const { return ev; }

    private:
        // Solver options
        Freccia::Options::DPR1EigenSolverOptions opt;
        
        // Input parameters
        int NS;
        int NR;

        DPR1Matrix<double> S; // Source DPR1 matrix
        
        DPR1Matrix<double> R; // Reduced DPR1 matrix
        Eigen::ArrayXd Rzsqr; // Square of R.z

        bool isRecasted = false; // Flag to indicate if R has been recasted to __float128
        DPR1Matrix<__float128> R_ld; // Reduced DPR1 matrix in __float128 datatype
        Eigen::Array<__float128, Eigen::Dynamic, 1> Rzsqr_ld; // __float128 version of Rzsqr

        // Sorting
        Permutation sort;

        // Type 1 deflation
        Freccia::Deflation::Type1<Freccia::Options::DPR1EigenSolverOptions> type1_deflation;
        
        // Type 2 deflation
        Freccia::Deflation::Type2 type2_deflation;
        
        // Output parameters
        Eigen::VectorXd ew; // Eigenwerte
        Eigen::MatrixXd ev; // Eigenvektoren

        // Functions
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

} // namespace DPR1
#endif