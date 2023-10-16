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

namespace Freccia::DPR1 {

class DPR1EigenSolver {
    public:
        // Constructor for ArrayXd
        DPR1EigenSolver(const Eigen::ArrayXd& D, 
                        const Eigen::ArrayXd& z,
                        double rho = 1.0,
                        const DPR1EigenSolverOptions& OPTIONS_IN = DPR1EigenSolverOptions()
                        ) 
        :
        sort(D.size()),
        type1_deflation(D.size()),
        type2_deflation(D.size()),
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

        // Return HH deflation matrix
        const HHDeflationMatrix & getHH() const { return HH; }
        
        // Return type 2 deflation indices
        const std::vector<unsigned int> & getType2DeflationIndices() { 
            return type2_deflation.nnzero();
        }

    private:
        // Solver options
        DPR1EigenSolverOptions opt;
        
        // Input parameters
        int NS;
        int NR;

        DPR1Matrix<double> S; // Source DPR1 matrix
        
        DPR1Matrix<double> R; // Reduced DPR1 matrix
        Eigen::ArrayXd Rzsqr; // Square of R.z

        bool isRecasted = false; // Flag to indicate if R has been recasted to __float128
        DPR1Matrix<__float128> R_ld; // Reduced DPR1 matrix in __float128 datatype
        Eigen::Array<__float128, Eigen::Dynamic, 1> Rzsqr_ld; // __float128 version of Rzsqr

        // Efficient storage of HH vectors for deflation
        HHDeflationMatrix HH;
        
        // Permutation and partitions
        Permutation sort;
        Partition type1_deflation;
        Partition type2_deflation;
        
        // Output parameters
        Eigen::VectorXd ew; // Eigenwerte
        Eigen::MatrixXd ev; // Eigenvektoren

        // Functions
        void eigh(const Eigen::ArrayXd& D, const Eigen::ArrayXd& z, double rho);
        void type1Deflation();
        void type2Deflation();
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