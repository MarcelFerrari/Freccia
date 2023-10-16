#ifndef FRECCIA_DPRK_DPRKEIGENSOLVER_H
#define FRECCIA_DPRK_DPRKEIGENSOLVER_H

// Eigen
#include <Eigen/Dense>

// Freccia
#include "Freccia/dpr1/DPR1EigenSolver.hpp"

namespace Freccia::DPRK {

    class DPRKEigenSolver {
        public:
        // Constructor for DPRK eigenproblem
        DPRKEigenSolver(const Eigen::Ref<const Eigen::VectorXd>& D, const Eigen::Ref<const Eigen::MatrixXd>& W, const double rho){
            DPRKeigh(D, W, rho);
        };

        // Constructor for rank K update of previously computed eigensystem
        DPRKEigenSolver(const Eigen::Ref<const Eigen::VectorXd>& D, const Eigen::Ref<const Eigen::MatrixXd>& W, const double rho, const Eigen::Ref<const Eigen::MatrixXd>& Q0){
            DPRKeigh(D, W, rho, Q0);
        };

        // Getters
        const Eigen::VectorXd& eigenvalues() const { return ew; };
        
        const Eigen::MatrixXd& eigenvectors() const { return ev; };
        
        private:
        // Member variables
        Eigen::VectorXd ew; // Eigenwerte
        Eigen::MatrixXd ev; // Eigenvektoren

        // Member functions
        void DPRKeigh(const Eigen::Ref<const Eigen::VectorXd>& D, const Eigen::Ref<const Eigen::MatrixXd>& W, const double rho);
        void DPRKeigh(const Eigen::Ref<const Eigen::VectorXd>& D, const Eigen::Ref<const Eigen::MatrixXd>& W, const double rho, const Eigen::Ref<const Eigen::MatrixXd>& Q0);
    };

} // End namespace Freccia::DPRK
#endif // FRECCIA_DPRK_DPRKEIGENSOLVER_H