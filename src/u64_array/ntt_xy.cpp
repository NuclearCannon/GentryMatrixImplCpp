#include "u64_array.hpp"
#include "ntt.hpp"


TwistedNtterXY64::TwistedNtterXY64(int n, u64 q, u64 qroot):
    n_(n), q_(q),
    zeta_pos_pows_(n),
    zeta_neg_pows_(n),
    std_ntter(n, q, qroot),
    mm(q)
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
    // std_ntter.ntt(dst.data(), buf.data());
    mm.batch_encode_inplace(buf);
    std_ntter.ntt_mont(dst.data(), buf.data());
    mm.batch_decode_inplace(dst);


}
void TwistedNtterXY64::intt(vec64& dst, const vec64& src) const
{
    vec64 buf(n_);
    // std_ntter.intt(buf.data(), src.data());
    vec64 src_encode = mm.batch_encode(src);
    std_ntter.intt_mont(buf.data(), src_encode.data());
    mm.batch_decode_inplace(buf);
    // let dst = buffer 逐位乘 zeta_neg_pows
    vec_mul(dst, buf, zeta_neg_pows_, q_);
}