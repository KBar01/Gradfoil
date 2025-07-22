#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"
#include "main_func.h"
#include <chrono>
#include <fstream>

using namespace std::chrono;



Real euc_norm(const Real* R, int size) {
    Real sum = 0.0;
    for (int i = 0; i < size; ++i) {
        sum += R[i] * R[i];
    }
    return std::sqrt(sum);
}


#ifdef USE_CODIPACK
bool solve_coupled(const Oper& oper, const Foil& foil, const Wake& wake,
    Param& param, Vsol& vsol, Isol& isol, Glob& glob, Trans&tdata, const bool force) {

    int nNewton = param.niglob;
    bool converged = false;
    constexpr int Rsize = 3*(Ncoords + Nwake);
    constexpr int Rallsize = 4*(Ncoords + Nwake);

    for (int i = 0; i < 1; ++i) {
        
        build_glob_RV(foil, vsol, isol, glob, param,tdata);

        solve_glob(foil, isol, glob, vsol, oper);

        update_state(oper, param, glob, vsol);
      
        clear_RV(glob, isol, vsol, foil, param);
        for (int entry = 0; entry < Rallsize; ++entry) {
            glob.R[entry] = 0;
        }

        stagpoint_move(isol, glob, foil, wake, vsol);
     
        update_transition(glob, vsol, isol, param, tdata, force);

        build_glob_RV(foil, vsol, isol, glob, param,tdata);

        Real residualNorm = euc_norm(glob.R, Rsize);

        if (residualNorm < param.rtol) {
            converged = true;
        }
    }

    return converged;
}

#else
bool solve_coupled(const Oper& oper, const Foil& foil, const Wake& wake,
    Param& param, Vsol& vsol, Isol& isol, Glob& glob, Trans&tdata, const bool force) {

    int nNewton = param.niglob;
    bool converged = false;
    constexpr int Rsize = 3*(Ncoords + Nwake);
    constexpr int Rallsize = 4*(Ncoords + Nwake);

    for (int i = 0; i < 50; ++i) {
        
        build_glob_RV(foil, vsol, isol, glob, param,tdata);
        
        Real residualNorm = euc_norm(glob.R, Rsize);
        
        if ((residualNorm < param.rtol) || (i == param.breakLoop)) {
            
            clear_RV(glob, isol, vsol, foil, param);
            converged = true;
            glob.convergenceIteration = i;
            break;
        }
        
        solve_glob(foil, isol, glob, vsol, oper);
        
        update_state(oper, param, glob, vsol);

        clear_RV(glob, isol, vsol, foil, param);
        for (int entry = 0; entry < Rallsize; ++entry) {
            glob.R[entry] = 0;
        }

        stagpoint_move(isol, glob, foil, wake, vsol);
        
        update_transition(glob, vsol, isol, param, tdata, force);

    }

    return converged;
}

#endif