#include "u64_array.hpp"
#include "ntt.hpp"


TwistedNtterXY64::TwistedNtterXY64(int n, u64 q, u64 qroot):
    n_(n), q_(q),
    zeta_pos_pows_(n),
    zeta_neg_pows_(n),
    std_ntter(n, q, qroot),
    mm(q),
    buf_n_(n)
{
    // 检查n：必须是power of 2
    assert(is_power_of_two(n));
    assert(q%(4*n)==1);
    u64 zeta = mod_pow(qroot, (q-1)/(4*n), q);
    // 检查zeta合法性
    u64 t = mod_pow(zeta, n*2, q);  // 理论上应该是-1
    assert(t+1==q);
    // zeta合法
    // 生成 zeta_pos_pows_
    get_powers(zeta_pos_pows_, zeta, n, q);
    // 生成 zeta_neg_pows_ninv_
    u64 zeta_inv = mod_inv(zeta, q);
    get_powers(zeta_neg_pows_, zeta_inv, n, q);

    zeta_pos_pows_mont_ = mm.batch_encode(zeta_pos_pows_);
    zeta_neg_pows_mont_ = mm.batch_encode(zeta_neg_pows_);

}

TwistedNtterXY64::~TwistedNtterXY64()
{
    // do nothing
}

void TwistedNtterXY64::ntt(vec64& dst, const vec64& src) const
{
    vec64 src_encode = mm.batch_encode(src);
    mm.vec_mul_mont(buf_n_, src_encode, zeta_pos_pows_mont_);
    std_ntter.ntt_mont(dst.data(), buf_n_.data());
    mm.batch_decode_inplace(dst);


}
void TwistedNtterXY64::intt(vec64& dst, const vec64& src) const
{
    vec64 src_encode = mm.batch_encode(src);
    std_ntter.intt_mont(buf_n_.data(), src_encode.data());
    // let dst = buffer 逐位乘 zeta_neg_pows
    mm.vec_mul_mont(dst, buf_n_, zeta_neg_pows_mont_);
    mm.batch_decode_inplace(dst);
}

void TwistedNtterXY64::ntt_mont(vec64& dst, const vec64& src) const
{
    mm.vec_mul_mont(buf_n_, src, zeta_pos_pows_mont_);
    std_ntter.ntt_mont(dst.data(), buf_n_.data());
}
void TwistedNtterXY64::intt_mont(vec64& dst, const vec64& src) const
{
    std_ntter.intt_mont(buf_n_.data(), src.data());
    // let dst = buffer 逐位乘 zeta_neg_pows
    mm.vec_mul_mont(dst, buf_n_, zeta_neg_pows_mont_);
}