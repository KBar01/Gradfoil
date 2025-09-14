#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

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


#ifndef USE_CODIPACK
void solve_sys(Glob& glob) {
    constexpr int Nsize = 4 * (Ncoords + Nwake);

    //auto start_total = high_resolution_clock::now();
    // Map the dense matrix from raw data
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
        A_eigen(glob.R_V, Nsize, Nsize);

    Eigen::Map<const Eigen::Matrix<Real, RVdimension, 1, Eigen::ColMajor>>
        rhs_eigen(glob.R, Nsize, 1);

    // Convert dense matrix to sparse matrix
    Eigen::SparseMatrix<Real> A_sparse(Nsize, Nsize);
    A_sparse = A_eigen.sparseView();  // Converts dense to sparse

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

    //auto t2 = high_resolution_clock::now();
    
    //auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - start_total);
    //std::cerr << "sparse compute: " << diff.count() << " " << std::endl;
}

#else

/*
template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename Type>
void func(Matrix<Type> const& A, Vector<Type> const& rhs, Vector<Type>& sol) {
    
    sol = A.colPivHouseholderQr().solve(rhs);
}
 
template<typename Number>
struct EigenSolver : public codi::EigenLinearSystem<Number, Matrix, Vector> {
  public:
 
    using Base = codi::EigenLinearSystem<Number, Matrix, Vector>;  
    using MatrixReal = typename Base::MatrixReal;                  
    using VectorReal = typename Base::VectorReal;                  
 
    void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x) {
        func(*A, *b, *x);
    }
};

void solve_sys(Glob&glob){
    
    
    constexpr int Nsize = 4 * (Ncoords + Nwake);
    

    Matrix<Real> A(Nsize, Nsize);
    Vector<Real> rhs(Nsize);
    Vector<Real> sol(Nsize);

    // Map the raw data to the Eigen matrices
    for (int i = 0; i < Nsize; ++i) {
        for (int j = 0; j < Nsize; ++j) {
            A(i, j) = glob.R_V[colMajorIndex(i, j, Nsize)];
        }
        rhs(i) = glob.R[i];
    }

    // Note: x here is a Eigen matrix 
    codi::solveLinearSystem(EigenSolver<Real>(), A, rhs, sol);
    
    for (int i=0;i<Nsize;++i){
        glob.dU[i] = -sol(i);
    }
}
*/

//#ifdef USE_CODIPACK

// === Drop-in sparse solver with full custom gradients for Codipack ===
// Place inside the #else branch where you currently have the dense QR solve.
// Requires Eigen (SparseLU) and Codipack headers already included.




#include <vector>
#include <iostream>
#include <cassert>

// keep your colMajorIndex helper
inline int colMajorIndex(int i, int j, int n) { return i + j * n; }

// The main function to replace the current codi branch solve_sys.
// It builds the sparse matrix from glob.R_V (reading numerical values from Real),
// factorizes, solves, computes the local Jacobian (A^{-1} and A^{-1} * rhs),
// and registers the per-output derivatives with Codipack.
void solve_sys(Glob &glob) {

    // sizes
    constexpr int Nsys = Ncoords + Nwake;
    const int Nsize = 4 * Nsys;

    // Helper to extract a plain double value from Real (works for both AD and non-AD Real).
    auto getValueDouble = [&](const Real &r) -> double {
    #ifdef USE_CODIPACK
        // In AD build Real is a Codipack type that exposes getValue()
        return r.getValue();
    #else
        return static_cast<double>(r);
    #endif
    };

    // 1) Build sparse matrix A from glob.R_V using triplets (only nonzeros)
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(4096); // adjust heuristic if you can estimate nnz

    for (int col = 0; col < Nsize; ++col) {
        for (int row = 0; row < Nsize; ++row) {
            int flat = colMajorIndex(row, col, Nsize);
            double val = getValueDouble(glob.R_V[flat]);
            if (val != 0.0) {
                triplets.emplace_back(row, col, val);
            }
        }
    }

    Eigen::SparseMatrix<double> A_sparse(Nsize, Nsize);
    A_sparse.setFromTriplets(triplets.begin(), triplets.end());

    // 2) Build rhs vector (double)
    Eigen::VectorXd rhs(Nsize);
    for (int i = 0; i < Nsize; ++i) rhs(i) = getValueDouble(glob.R[i]);

    // 3) Factorize once with SparseLU (reuse for multiple solves)
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A_sparse);
    solver.factorize(A_sparse);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Sparse factorization failed (forward)!\n";
        // You may want to set glob.dU to zero or handle error gracefully.
        for (int i = 0; i < Nsize; ++i) glob.dU[i] = (Real)0.0;
        return;
    }

    // 4) Solve y = A^{-1} * rhs
    Eigen::VectorXd y = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Sparse solve failed (forward)!\n";
        for (int i = 0; i < Nsize; ++i) glob.dU[i] = (Real)0.0;
        return;
    }

    // 5) Compute solution x = -y and write into glob.dU (as Real)
    Eigen::VectorXd x = -y;
    for (int i = 0; i < Nsize; ++i) {
        glob.dU[i] = (Real)x(i); // constructs AD variable when in AD build
    }

    // ------------------------------
    // 6) Compute A^{-1} to build local Jacobian entries:
    //    we need Ainv = A^{-1} so we can compute:
    //      ∂x_i/∂b_j = - (A^{-1})_{i,j}
    //      ∂x_i/∂A_{p,q} = (A^{-1})_{i,p} * y_q
    //    We compute Ainv by solving A * E = I (multiple RHS) using the same factorization.
    // ------------------------------

    // Build identity matrix as dense double (Nsize x Nsize) and solve in one call
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(Nsize, Nsize);
    Eigen::MatrixXd Ainv = solver.solve(I);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solve for Ainv failed!\n";
        // fallback: we could compute columns one-by-one, but we stop here for clarity
        for (int i = 0; i < Nsize; ++i) {
            glob.dU[i] = (Real)0.0;
        }
        return;
    }

    // ------------------------------
    // 7) Register per-output local derivatives with Codipack using StatementPushHelper.
    //    For each output i (glob.dU[i]) we push:
    //      - each b_j with derivative = ∂x_i/∂b_j = - Ainv(i,j)
    //      - each *nonzero* A_{p,q} (matching glob.R_V layout) with derivative = Ainv(i,p) * y(q)
    //    This gives Codipack exactly the Jacobian needed at this linear solve node.
    // ------------------------------

    // To speed up the inner loops we cache the nonzero triplet info we already extracted.
    // Build arrays of nz rows/cols/flat indices for fast iteration
    std::vector<int> nz_rows; nz_rows.reserve(triplets.size());
    std::vector<int> nz_cols; nz_cols.reserve(triplets.size());
    std::vector<int> nz_flat; nz_flat.reserve(triplets.size());
    for (const auto &t : triplets) {
        nz_rows.push_back((int)t.row());
        nz_cols.push_back((int)t.col());
        nz_flat.push_back(colMajorIndex((int)t.row(), (int)t.col(), Nsize));
    }

    // For each output i create a statement and push all input arguments with the correct scalar partials
    // IMPORTANT: StatementPushHelper must use the Codipack reverse real type
    for (int i = 0; i < Nsize; ++i) {

        // Create a push helper for this output
        codi::StatementPushHelper<codi::RealReverse> ph;
        ph.startPushStatement();

        // 7a) push each RHS element: glob.R[j] with derivative -Ainv(i,j)
        for (int j = 0; j < Nsize; ++j) {
            double deriv_b = - Ainv(i, j);          // ∂x_i / ∂b_j
            ph.pushArgument(glob.R[j], deriv_b);
        }

        // 7b) push each nonzero A element from glob.R_V:
        // derivative = ∂x_i / ∂A_{p,q} = Ainv(i,p) * y(q)
        for (size_t k = 0; k < nz_rows.size(); ++k) {
            int p = nz_rows[k];   // row index
            int q = nz_cols[k];   // col index
            int flat = nz_flat[k]; // flattened index into glob.R_V
            double deriv_A = Ainv(i, p) * y(q);    // Ainv(i,p) * y_q
            ph.pushArgument(glob.R_V[flat], deriv_A);
        }

        // 7c) finish the statement for output glob.dU[i]
        // supply the current (double) value of x(i) so Codipack records it.
        ph.endPushStatement(glob.dU[i], x(i));
    }

    // Done. The AD tape now has a single "atomic" statement per output x_i where the local Jacobian
    // entries are exact numeric values computed above. Codipack will use those numbers during reverse,
    // so we avoided tracing the dense factorization internals.
}

#endif

void writeArrayToCSV(const std::string& filename, const double* array, int size) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    for (int i = 0; i < size; ++i) {
        file << array[i] << "\n";
    }

    file.close();
}


void solve_glob(const Foil&foil, const Isol&isol, Glob& glob, Vsol& vsol, const Oper& oper) {
    
    
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
            glob.R_V[colMajorIndex(rowStart+row,colindex,4*Nsys)] = (row == col ? 1.0 : 0.0) - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ds[col];
        }
    }

    //all disp thickness indices
    for (int col=0;col<Nsys;++col){

        int colindex = 4*col + 1;
        for (int row = 0;row<Nsys;++row){
            glob.R_V[colMajorIndex(rowStart+row,colindex,4*Nsys)] =  - vsol.ue_m[colMajorIndex(row,col,Nsys)]*ue[col];
        }
    }

    solve_sys(glob);
}


