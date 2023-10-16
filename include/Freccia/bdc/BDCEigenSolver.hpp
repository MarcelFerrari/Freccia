#include <Eigen/Dense>
#include <vector>

#ifndef TBD_H
#define TBD_H

namespace Freccia::Banded {
class BDCEigenSolver {
    using SVDTuple = std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd>;
    using EigenTuple = std::tuple<unsigned int, unsigned int>;
    
    public:
        // Input matrix A is a banded matrix with bandwidth f
        BDCEigenSolver(const Eigen::MatrixXd &A, const unsigned int f_in){
            block_divide_and_conquer(A, f_in); // Call block divide and conquer algorithm
        };

        // Getters
        Eigen::VectorXd eigenvalues() const { return D;}
        Eigen::MatrixXd eigenvectors() const { return Q;}
        
    private:
        // Iput variables
        unsigned int N; // Size of A
        unsigned int f; // Size of off-diagonal blocks
        unsigned int p; // Number of diagonal blocks

        // Eigenvalues and eigenvectors
        std::vector<SVDTuple> svdBlocks;
        Eigen::VectorXd D;
        Eigen::MatrixXd Q;

        // Functions
        void block_divide_and_conquer(const Eigen::MatrixXd &A, const unsigned int f_in);
        EigenTuple mergeBlocks(EigenTuple& A1, EigenTuple& A2);
        
   };
}

#endif
