#ifndef FRECCIA_OPTIONS_H
#define FRECCIA_OPTIONS_H

#include "Freccia/utils/utils.hpp"

namespace Freccia::Options {
    // Struct for DPR1 Eigensolver options
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

    // Arrowhead solver uses similar options as DPR1, but different environment variables
    struct ArrowheadEigenSolverOptions {
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
        ArrowheadEigenSolverOptions() {}
        
        // Copy constructor
        ArrowheadEigenSolverOptions(const DPR1EigenSolverOptions& other)
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
                read_env("FRECCIA_ARROW_KZ_TOL", KZ_TOL);
                read_env("FRECCIA_ARROW_KB_TOL", KB_TOL);
                read_env("FRECCIA_ARROW_KNU_TOL", KNU_TOL);
                read_env("FRECCIA_ARROW_RHO_TOL", RHO_TOL);
                read_env("FRECCIA_ARROW_K_LAMBDA_TOL", K_LAMBDA_TOL);
                read_env("FRECCIA_ARROW_BISECT_MAX_ITER", BISECT_MAX_ITER);
                read_env("FRECCIA_ARROW_RECOMPUTE_MAX_ITER", RECOMPUTE_MAX_ITER);
                read_env("FRECCIA_ARROW_ABS_ZERO_TOL", ABS_ZERO_TOL);
                read_env("FRECCIA_ARROW_REL_ZERO_TOL", REL_ZERO_TOL);
            }
        }
    };
}

#endif