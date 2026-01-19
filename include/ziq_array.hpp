#pragma once
#include <vector>
#include <stdexcept>
#include "flints.hpp"
#include "twisted_ntter.hpp"


class MatmulContext {
    int n_;
    fmpz_t q_;
    fmpz_mod_mat_t A_mat, B_mat, C_mat; // 矩阵乘法时使用的缓冲区
public:
    MatmulContext(int n, fmpz_t q);
    ~MatmulContext();

    // C = A @ B.T  (mod q)
    void matmul_transpose(fmpz* C, const fmpz* A, const fmpz* B);
};


class ZiqArrayContext {
private:
    int n_, p_;
    int nn_, pnn_, size_;
    fmpz_t q_, I_, I_inv_;
    fmpz_mod_ctx_t q_ctx_;
    TwistedNtterXY *ntter_p, *ntter_n;    // 这两个ntter会分别负责模X^n-I的和模X^n+I的
    TwistedNtterW *ntter_w;
    MatmulContext *mm_ctx;

    fmpz_vector *buf_p, *buf_n, *buf_half, *buf_size;
public:

    ZiqArrayContext(int n, int p, int g, fmpz_t q, fmpz_t zeta, fmpz_t eta);

    ~ZiqArrayContext();

    fmpz_vector iw_ntt(const fmpz_vector& src) const;
    fmpz_vector iw_intt(const fmpz_vector& src) const;

    fmpz_vector xy_ntt(const fmpz_vector& src) const;
    fmpz_vector xy_intt(const fmpz_vector& src) const;

    fmpz_vector circledast(const fmpz_vector& A, const fmpz_vector& B) const;

};

