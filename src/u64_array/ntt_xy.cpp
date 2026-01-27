#include "u64_array.hpp"
#include "ntt.hpp"


TwistedNtterXY64::TwistedNtterXY64(int n, u64 q, u64 zeta):
    n_(n), q_(q),
    zeta_pos_pows_(n),
    zeta_neg_pows_ninv_(n),
    omega_pos_pows_(n),
    omega_neg_pows_(n)

{
    // 检查n：必须是power of 2
    assert(is_power_of_two(n));
    // 检查zeta合法性
    u64 t = mod_pow(zeta, n*2, q);  // 理论上应该是-1
    assert(t+1==q);
    // zeta合法
    // 生成 zeta_pos_pows_
    get_powers(zeta_pos_pows_, zeta, n, q);
    // 生成 zeta_neg_pows_ninv_
    u64 ninv = mod_inv(n, q);
    u64 zeta_inv = mod_inv(zeta, q);
    get_powers(zeta_neg_pows_ninv_, zeta_inv, n, q);
    vec_scalar_mul(zeta_neg_pows_ninv_, zeta_neg_pows_ninv_, ninv, q);
    // 生成 omega_pos_pows_
    u64 omega = mod_pow(zeta, 4, q);
    get_powers(omega_pos_pows_, omega, n, q);
    u64 omega_inv = mod_inv(omega, q);
    get_powers(omega_neg_pows_, omega_inv, n, q);
}

TwistedNtterXY64::~TwistedNtterXY64()
{
    // do nothing
}

void TwistedNtterXY64::ntt(vec64& dst, const vec64& src) const
{
    vec64 buf(n_);
    // let buf = src 逐位乘 zeta_pos_pows
    vec_mul(buf, src, zeta_pos_pows_, q_);
    ntt_standard_64_with_roots(buf.data(), dst.data(), omega_pos_pows_.data(), n_, q_);


}
void TwistedNtterXY64::intt(vec64& dst, const vec64& src) const
{
    vec64 buf(n_);
    // let buffer = NTT(buffer, omega^{-1}, n, ctx)
    ntt_standard_64_with_roots(src.data(), buf.data(), omega_neg_pows_.data(), n_, q_);
    // let dst = buffer 逐位乘 zeta_neg_pows_ninv
    vec_mul(dst, buf, zeta_neg_pows_ninv_, q_);
}