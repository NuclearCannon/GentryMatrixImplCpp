#pragma once
#include <vector>
#include <stdexcept>
#include "flints.hpp"
#include "twisted_ntter.hpp"

class ZiqArrayContext {
private:
    int n_, p_;
    int nn_, pnn_, size_;
    fmpz_t q_, I_, I_inv_;
    fmpz_mod_ctx_t q_ctx_;
    TwistedNtterXY *ntter_p, *ntter_n;    // 这两个ntter会分别负责模X^n-I的和模X^n+I的
    TwistedNtterW *ntter_w;

    fmpz_vector *buf_p, *buf_n, *buf_half, *buf_size;
public:

    ZiqArrayContext(int n, int p, int g, fmpz_t q, fmpz_t zeta, fmpz_t eta);

    ~ZiqArrayContext();

    fmpz_vector iw_ntt(const fmpz_vector& src) const;
    fmpz_vector iw_intt(const fmpz_vector& src) const;

    fmpz_vector xy_ntt(const fmpz_vector& src) const;
    fmpz_vector xy_intt(const fmpz_vector& src) const;



};

