#ifndef FRECCIA_TYPE_1_DEFLATION_H
#define FRECCIA_TYPE_1_DEFLATION_H

#include <Eigen/Dense>
#include "Freccia/matrix/matrix.hpp"
#include "Freccia/options/options.hpp"
#include "Freccia/utils/numeric.hpp"

// The deflation procedures are templated to work both for DPR1 and Arrowhead solvers
namespace Freccia::Deflation {
    class Type1{
        public:
            // Constructor
            Type1() {}

            // Deflate interface
            void deflate(const Eigen::ArrayXd& D, Eigen::ArrayXd& z, const double ABS_ZERO_TOL_IN, const double REL_ZERO_TOL_IN){
                ABS_ZERO_TOL = ABS_ZERO_TOL_IN;
                REL_ZERO_TOL = REL_ZERO_TOL_IN;
                compute_type1_deflation(D, z); // Call type 1 deflation
            };

            // Return HH deflation matrix
            const HHDeflationMatrix & getHH() const { return HH; }

            // Return permutation
            const Partition & getPartition() const { return type1_deflation; }
        
        private:
            // Need to store partitioning and HH matrix for type 1 deflation
            Partition type1_deflation;
            HHDeflationMatrix HH;
            double ABS_ZERO_TOL;
            double REL_ZERO_TOL;
            
            void compute_type1_deflation(const Eigen::ArrayXd& D_in, Eigen::ArrayXd& z_in){
                // Partition values in z_in
                type1_deflation.partition(z_in, ABS_ZERO_TOL);

                // Apply type 1 deflation
                const std::vector<unsigned int>& nnzero = type1_deflation.nnzero();
                
                // Store eigen IndexedView to modify z in place
                const auto D = D_in(nnzero);
                auto z = z_in(nnzero); 

                // Compute multiplicity of D matrix
                std::vector<std::pair<unsigned int, unsigned int>> multiplicity;

                // Start a scope to localize the variables used in this block
                {
                    // Here N refers to the size of the reduced problem
                    unsigned int N = D.size();
                    unsigned int i = 0;

                    // Iterate over D
                    while(i < N){
                        unsigned int j = i + 1;
                        // Extract the current eigenvalue
                        double lam = D(i);
                        
                        // Increment 'j' while subsequent elements are approximately equal to the current eigenvalue.
                        // This will effectively group together and count elements of the same value (up to a small tolerance).
                        while(j < N && (Freccia::Numeric::isclose(D(j), lam, ABS_ZERO_TOL, REL_ZERO_TOL))){
                            ++j;
                        } 

                        // Store starting index and multiplicity of the current eigenvalue if it is degenerate
                        if(j - i > 1){ //&& std::abs(S.D(i)) > 1e-16){ // Deflate only if the eigenvalue is non-zero
                            multiplicity.push_back(std::make_pair(i, j - i));
                        }

                        // Update 'i' to continue with the next distinct eigenvalue
                        i = j;
                    }
                }

                // Iterate over the multiplicity vector which contains pairs of (start_index, multiplicity) = (i, m)
                for(const auto & [i, m]: multiplicity){
                    // Calculate the norm of the current z segment
                    double znorm = z.segment(i, m).matrix().norm();

                    // Extract the segment from z_sorted corresponding to this eigenvalue group
                    Eigen::VectorXd u = z.segment(i, m).matrix();

                    // Add the norm of the segment to the first element of 'u' to create the reflection vector
                    // if we choose u = z + |z| * e_0, then H*z = -|z| * e_0,
                    // alternatively, we can choose u = z - |z| * e_0 and obtain H*z = |z| * e_0
                    u(0) += znorm;
                    u.normalize(); // Normalize the reflection vector
            
                    // Push block to Householder matrix using the reflection vector 'u'
                    // The HHDeflationMatrix HH is defined in the class definition!
                    HH.push(i, m, std::move(u)); // Pass an rvalue reference to 'u' to avoid copying

                    // Deflate z by implicitly applying the Householder reflection: 
                    // the first element becomes the norm, and the rest of the elements become zero.
                    z(i) = -znorm; // The sign is important here! 
                    z.segment(i + 1, m - 1).setZero();
                }
            }
    };
}

#endif // FRECCIA_TYPE_1_DEFLATION_H