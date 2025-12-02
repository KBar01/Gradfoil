#include <iostream>
#include <cmath>
#include <Eigen/Dense>


#include <Eigen/LU>
#include "real_type.h"
#include "panel_funcs.h"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"


#include <vector>
#include <cmath>
#include <algorithm> // for std::max
#include <iostream>

#include <chrono>
#include <fstream>

using namespace std::chrono;




// === Drop-in sparse solver with full custom gradients for Codipack ===
// Place inside the #else branch where you currently have the dense QR solve.
// Requires Eigen (SparseLU) and Codipack headers already included.

#include <vector>
#include <iostream>
#include <cassert>

#include <vector>
#include <iostream>
#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>






# ifdef USE_CODIPACK
// helper for column-major indexing
inline int colMajorIndex(int i, int j, int n) { return i + j * n; }

// 1. Type aliases for sparse Eigen
template<typename T>
using MatrixSparse = Eigen::SparseMatrix<T>;
template<typename T>
using Vector      = Eigen::Matrix<T, Eigen::Dynamic, 1>;

// 2. Your own solver for the numeric step
template<typename T>
void sparseSolveFunc(MatrixSparse<T> const& A, Vector<T> const& rhs, Vector<T>& sol) {
    // choose your factorization; SparseLU works for general unsymmetric
    Eigen::SparseLU<MatrixSparse<T>> lu;
    lu.compute(A);
    sol = lu.solve(rhs);
}

// 3. Wrap in CoDiPack's sparse linear system
template<typename Number>
struct SparseEigenSolver
  : public codi::SparseEigenLinearSystem<Number, MatrixSparse, Vector> {

    using Base       = codi::SparseEigenLinearSystem<Number, MatrixSparse, Vector>;
    using MatrixReal = typename Base::MatrixReal;  // numeric (Real) matrix
    using VectorReal = typename Base::VectorReal;  // numeric (Real) vector

    void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x) {
        sparseSolveFunc(*A, *b, *x);  // just delegate to your numeric routine
    }
};

// 4. Your driver
void solve_sys(Glob &glob) {
    constexpr int Nsize = 4 * (Ncoords + Nwake);

    // build sparse matrix from your glob arrays
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(glob.R_V_latest);
    for (int k = 0; k < glob.R_V_latest; ++k) {
        triplets.emplace_back(glob.R_V_rows[k],
                              glob.R_V_cols[k],
                              glob.R_V_vals[k]);
    }

    MatrixSparse<Real> A(Nsize, Nsize);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Vector<Real> rhs(Nsize);
    for (int i = 0; i < Nsize; ++i) rhs(i) = glob.R[i];

    Vector<Real> sol(Nsize);

    // 5. Let CoDiPack handle AD by calling its wrapper:
    codi::solveLinearSystem(SparseEigenSolver<Real>(), A, rhs, sol);

    // 6. Write solution back into your state
    for (int i = 0; i < Nsize; ++i) {
        glob.dU[i] = -sol(i);
    }
}

#else

void solve_sys(Glob& glob) {
    constexpr int Nsize = 4 * (Ncoords + Nwake);

    // === 1. Build triplets from glob arrays ===
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(glob.R_V_latest);
    for (int k = 0; k < glob.R_V_latest; ++k) {
        int row = glob.R_V_rows[k];
        int col = glob.R_V_cols[k];
        Real val = glob.R_V_vals[k];
        if (val != Real(0)) {
            triplets.emplace_back(row, col, val);
        }
    }

    // === 2. Fill sparse matrix from triplets ===
    Eigen::SparseMatrix<Real> A_sparse(Nsize, Nsize);
    A_sparse.setFromTriplets(triplets.begin(), triplets.end());
    
    Eigen::Map<const Eigen::Matrix<Real, RVdimension, 1, Eigen::ColMajor>>
    rhs_eigen(glob.R, Nsize, 1);

    // Use SparseLU solver
    Eigen::SparseLU<Eigen::SparseMatrix<Real>> sparse_solver;
    sparse_solver.compute(A_sparse);

    if(sparse_solver.info() != Eigen::Success) {
        std::cerr << "Sparse solver failed during factorization!\n";
        return; // or handle the error appropriately
    }
    
    Eigen::Matrix<Real, RVdimension, 1> x = -sparse_solver.solve(rhs_eigen);

    // Map the solution to output vector
    Eigen::Map<Eigen::Matrix<Real, RVdimension, 1, Eigen::ColMajor>>
        x_eigen(glob.dU, Nsize, 1);
    x_eigen = x;
}
#endif


void solve_glob(const Foil&foil, const Isol&isol, Glob& glob, Vsol& vsol, const Oper& oper, const int doSolve) {
    
    
    constexpr int Nsys = Ncoords+Nwake;

    // Step 1: Modify ue array to avoid 0 or negative
    int nrows = 4; // Since U is shaped (4, Nsys) in column-major
    Real ue[Nsys] = {0};
    Real uemax = 0.0;
    for (int i = 0; i < Nsys; ++i){
        uemax = std::max(uemax, std::abs(glob.U[colMajorIndex(3,i,4)]));
    }
    for (int i = 0; i < Nsys; ++i){
        ue[i] = std::max(glob.U[colMajorIndex(3,i,4)], 1e-10*uemax);
    }

    // Step 2: Get ueinv
    Real ueinv[Nsys]={0};
    get_ueinv(isol,ueinv);

    // Step 3: Build Residual R
    Real ds[Nsys];
    for (int i = 0; i < Nsys; ++i) {ds[i] = glob.U[colMajorIndex(1, i, 4)];}

    Real tempRHS[Nsys];
    cnp::mul<Nsys>(ds,ue,tempRHS); // ds*ue

    Real* Rpointer = &glob.R[3*Nsys] ; 
    cnp::matmat_mul<Nsys,Nsys,1>(vsol.ue_m,tempRHS,Rpointer);

    for (int i = 0; i < Nsys; ++i){
        Rpointer[i] = ue[i] - (ueinv[i] + Rpointer[i]);
    }

    // all edge velocity indices
    int rowStart = 3*Nsys;
    for (int col=0;col<Nsys;++col){

        int colindex = 4*col + 3;
        for (int row = 0;row<Nsys;++row){
            
            //glob.R_V[colMajorIndex(rowStart+row,colindex,4*Nsys)] = (row == col ? 1.0 : 0.0) - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ds[col];
            
            Real zero = 0.0;
            glob.R_V_vals[glob.R_V_latest] = (row == col ? 1.0 : zero) - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ds[col];
            glob.R_V_rows[glob.R_V_latest] = rowStart+row;
            glob.R_V_cols[glob.R_V_latest] = colindex;
            glob.R_V_latest += 1 ;
        }
    }

    //all disp thickness indices
    for (int col=0;col<Nsys;++col){

        int colindex = 4*col + 1;
        for (int row = 0;row<Nsys;++row){
            //glob.R_V[colMajorIndex(rowStart+row,colindex,4*Nsys)] =  - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ue[col];
            glob.R_V_vals[glob.R_V_latest] = - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ue[col];
            glob.R_V_rows[glob.R_V_latest] = rowStart+row;
            glob.R_V_cols[glob.R_V_latest] = colindex;
            glob.R_V_latest += 1 ;
        }
    }

    if (doSolve) {
        solve_sys(glob);
    }
}


