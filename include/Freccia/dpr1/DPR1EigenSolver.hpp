#ifndef FRECCIA_DPR1EIGENSOLVER_H
#define FRECCIA_DPR1EIGENSOLVER_H

// Standard lib
#include <tuple>

// Eigen
#include <Eigen/Dense>

// Freccia
#include "Freccia/matrix/matrix.hpp"
#include "Freccia/utils/utils.hpp"

namespace Freccia::DPR1 {

// Struct for Eigensolver options
struct DPR1EigenSolverOptions {
    double rho = 1.0;
    double ABS_ZERO_TOL = 1e-13;
    double REL_ZERO_TOL = 1e-13;
    unsigned int KZ_TOL = 1e3;
    double KB_TOL = 1e3;
    double KNU_TOL = 1e3;
    double RHO_TOL = 1e3;
    double K_LAMBDA_TOL = 1e3;
    unsigned int BISECT_MAX_ITER = 100;
    unsigned int RECOMPUTE_MAX_ITER = 3;
    bool CHECK_ENV = true;
    bool RECONSTRUCT_Q = true;

    // Default constructor
    DPR1EigenSolverOptions() {}
    
    // Copy constructor
    DPR1EigenSolverOptions(const DPR1EigenSolverOptions& other)
        : rho(other.rho), ABS_ZERO_TOL(other.ABS_ZERO_TOL),
          REL_ZERO_TOL(other.REL_ZERO_TOL), KZ_TOL(other.KZ_TOL),
          KB_TOL(other.KB_TOL), KNU_TOL(other.KNU_TOL),
          RHO_TOL(other.RHO_TOL), K_LAMBDA_TOL(other.K_LAMBDA_TOL),
          BISECT_MAX_ITER(other.BISECT_MAX_ITER), 
          RECOMPUTE_MAX_ITER(other.RECOMPUTE_MAX_ITER),
          CHECK_ENV(other.CHECK_ENV), RECONSTRUCT_Q(other.RECONSTRUCT_Q)
    {}

    void loadEnv() {
        if(CHECK_ENV){
            read_env("FRECCIA_DPR1_KZ_TOL", KZ_TOL);
            read_env("FRECCIA_DPR1_KB_TOL", KB_TOL);
            read_env("FRECCIA_DPR1_KNU_TOL", KNU_TOL);
            read_env("FRECCIA_DPR1_RHO_TOL", RHO_TOL);
            read_env("FRECCIA_DPR1_K_LAMBDA_TOL", K_LAMBDA_TOL);
            read_env("FRECCIA_DPR1_BISECT_MAX_ITER", BISECT_MAX_ITER);
            read_env("FRECCIA_DPR1_RECOMPUTE_MAX_ITER", RECOMPUTE_MAX_ITER);
            read_env("FRECCIA_DPR1_ABS_ZERO_TOL", ABS_ZERO_TOL);
            read_env("FRECCIA_DPR1_REL_ZERO_TOL", REL_ZERO_TOL);
        }
    }
};

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