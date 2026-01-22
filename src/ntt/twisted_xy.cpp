#include "twisted_ntter.hpp"
#include "ntt.hpp"


TwistedNtterXY::TwistedNtterXY(int n, const fmpz_t q, const fmpz_t zeta):
    n_(n),
    zeta_pos_pows_(n), 
    zeta_neg_pows_ninv_(n),
    buffer_(n)
{
    assert(n > 0);
    // 初始化各个数字
    fmpz_init_set(q_, q);
    fmpz_mod_ctx_init(q_ctx_, q_); // 预计算模数 q 的优化参数

    fmpz_init(I_);
    fmpz_init_set(zeta_, zeta);
    fmpz_init(omega_);
    fmpz_init(omega_inv_);
    // I=zeta^n
    fmpz_mod_pow_ui(I_, zeta, n, q_ctx_);
    fmpz_mod_pow_ui(omega_, zeta, 4, q_ctx_);
    fmpz_invmod(omega_inv_, omega_, q_);

    // 生成zeta_pos_pows_
    fmpz_set_ui(zeta_pos_pows_[0], 1);
    for(int i=1;i<n;i++)
    {
        fmpz_mod_mul(zeta_pos_pows_[i], zeta_pos_pows_[i-1], zeta_, q_ctx_);
    }
    // zeta_neg_pows_ninv_
    fmpz_t n_mpz, n_inv, zeta_inv;
    fmpz_init_set_ui(n_mpz, n);
    fmpz_init(n_inv);
    fmpz_init(zeta_inv);

    int r = fmpz_invmod(n_inv, n_mpz, q_);
    assert(r);

    r = fmpz_invmod(zeta_inv, zeta_, q_);
    assert(r);

    fmpz_set(zeta_neg_pows_ninv_[0], n_inv);
    for(int i=1;i<n;i++)
    {
        fmpz_mod_mul(zeta_neg_pows_ninv_[i], zeta_neg_pows_ninv_[i-1], zeta_inv, q_ctx_);
    }

    fmpz_clear(n_mpz);
    fmpz_clear(n_inv);
    fmpz_clear(zeta_inv);
}

TwistedNtterXY::~TwistedNtterXY()
{
    fmpz_clear(q_);
    fmpz_clear(I_);
    fmpz_clear(zeta_);
    fmpz_clear(omega_);
    fmpz_clear(omega_inv_);
    fmpz_mod_ctx_clear(q_ctx_);
    // fmpz_vector自己会析构，不用管它

}
void TwistedNtterXY::ntt(const fmpz_vector& src, fmpz_vector& dst) 
{
    // let buffer = src 逐位乘 zeta_pos_pows
    _fmpz_mod_vec_mul(buffer_.raw(), src.raw(), zeta_pos_pows_.raw(), n_, q_ctx_);
    // let dst = NTT(buffer, omega, n, ctx)
    ntt_standard_flint(buffer_, dst, omega_, n_, q_ctx_);
    // 结束
}

void TwistedNtterXY::intt(const fmpz_vector& src, fmpz_vector& dst)
{
    // let buffer = NTT(buffer, omega^{-1}, n, ctx)
    ntt_standard_flint(src, buffer_, omega_inv_, n_, q_ctx_);
    // let dst = buffer 逐位乘 zeta_neg_pows_ninv
    _fmpz_mod_vec_mul(dst.raw(), buffer_.raw(), zeta_neg_pows_ninv_.raw(), n_, q_ctx_);
    // 结束
}
