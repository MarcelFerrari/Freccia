#ifndef FRECCIA_TYPE_2_DEFLATION_H
#define FRECCIA_TYPE_2_DEFLATION_H

#include <Eigen/Dense>
#include "Freccia/matrix/matrix.hpp"
#include "Freccia/options/options.hpp"

namespace Freccia::Deflation {
    class Type2 {
            public:
            // Constructor
            Type2() {}

            // Deflate interface
            void deflate(const Eigen::ArrayXd& D, const Eigen::ArrayXd& z, Eigen::VectorXd& ew, Eigen::MatrixXd& ev){
                compute_type2_deflation(D, z, ew, ev); // Call type 2 deflation
            };

            // Return partition
            const Partition & getPartition() const { return type2_deflation; }
        
        private:
            // Need to store partitioning and HH matrix for type 1 deflation
            Partition type2_deflation;
            
            void compute_type2_deflation(const Eigen::ArrayXd& D, const Eigen::ArrayXd& z, Eigen::VectorXd& ew, Eigen::MatrixXd& ev){
                // Partition zero and non zero values in S.z
                type2_deflation.partition(z);
                
                // Compute the deflated solution. For each index 'k' corresponding to a zero 
                // element in 'z', set the (k,k)-th element of matrix Q to 1 and the k-th 
                // element of vector L to the corresponding element in D_sorted. This effectively 
                // sets the diagonal elements of Q and L for zero elements, and leaves non-zero 
                // elements to be handled separately.
                for(unsigned int k: type2_deflation.zero()){
                    ev(k, k) = 1.;
                    ew(k) = D(k);
                }
            }
    };
}

#endif // FRECCIA_TYPE_1_DEFLATION_H