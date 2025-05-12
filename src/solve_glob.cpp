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


