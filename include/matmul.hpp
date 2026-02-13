#pragma once
#include "flints.hpp"
#include <cstdint>
#include "montgomery.hpp"


class MatmulContext {
    int n_;
    fmpz_t q_;
    mutable fmpz_mod_mat_t A_mat, B_mat, C_mat;
public:
    MatmulContext(int n, const fmpz_t q);
    ~MatmulContext();

    // C = A @ B.T  (mod q)
    void matmul_transpose(fmpz* C, const fmpz* A, const fmpz* B) const;

    void matmul_transpose_u64(uint64_t* C, const uint64_t* A, const uint64_t* B) const;

    fmpz_vector circledast_fmpz(const fmpz_vector& A, const fmpz_vector& B, size_t n, size_t p);

    void circledast_u64(uint64_t* dst, const uint64_t* A, const uint64_t* B, size_t n, size_t p);
};



void circledast_u64_gpu(uint64_t* dst, const uint64_t* A, const uint64_t* B, size_t n, size_t p, const MontgomeryMultiplier& mm);