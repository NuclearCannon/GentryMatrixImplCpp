#pragma once
#include "flints.hpp"
#include "uint64.hpp"

class MatmulContext {
    int n_;
    fmpz_t q_;
    mutable fmpz_mod_mat_t A_mat, B_mat, C_mat;
public:
    MatmulContext(int n, const fmpz_t q);
    ~MatmulContext();

    // C = A @ B.T  (mod q)
    void matmul_transpose(fmpz* C, const fmpz* A, const fmpz* B) const;
    fmpz_vector circledast(const fmpz_vector& A, const fmpz_vector& B, size_t n, size_t p);
};

fmpz_vector circledast(const fmpz_vector& A, const fmpz_vector& B, size_t n, size_t p, const MatmulContext& mc);
