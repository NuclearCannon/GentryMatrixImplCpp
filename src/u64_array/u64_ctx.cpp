#include "u64_array.hpp"
#include <cstring>
#include "math_utils.hpp"


U64Context::U64Context(int n, int p, uint64_t q, uint64_t root_q):
    n_(n),
    p_(p),
    nn_(n*n),
    pnn_((p-1)*n*n),
    size_(2*(p-1)*n*n),
    q_(q),
    mm_(q),
    buf_size_(2*(p-1)*n*n),
    bufn_(n)
{
    // 计算zeta. zeta是4n阶本原单位根
    assert((q-1)%((uint64_t)(4*n)) == 0);;
    // 计算eta. eta是p阶本原单位根
    assert((q-1)%((uint64_t)p) == 0);
    // 生成各ntter对象
    ntter_p = std::make_unique<TwistedNtterXY64>(n,q,root_q);
    ntter_n = std::make_unique<TwistedNtterXY64>(n,q,mod_inv(root_q, q));
    ntter_w = std::make_unique<TwistedNtterW64>(p,q,root_q);
    // 计算I及其逆元
    // I = zeta^n
    // zeta = qroot ^ ((q-1)/4n)
    // => I=qroot^((q-1)/4)
    I_ = mod_pow(root_q, (q-1)/4, q);
    I_inv_ = mod_pow(I_, 3, q);
    assert(mod_mul(I_, I_inv_, q) == 1);

    I_mont_ = mm_.encode(I_);
    I_inv_mont_ = mm_.encode(I_inv_);

    inv2_mont = mm_.encode(mod_inv(2, q));

}

U64Context::~U64Context()
{
    // do nothing
}

void U64Context::iw_ntt(vec64& dst, const vec64& src) const
{
    // ntt-I
    for(int i=0, j=pnn_;i<pnn_;i++,j++)
    {
        uint64_t real_i = src[i];
        uint64_t image_I_i = mm_.mul(src[j], I_mont_);
        dst[i] = mod_add(real_i, image_I_i, q_);
        dst[j] = mod_sub(real_i, image_I_i, q_);
    }
    // 现在dst是I-ntted

    // 转置：(2, p-1, N^2) -> (N^2, 2, p-1)
    vec64& buf = buf_size_;
    transpose_rect_restrict(buf.data(), dst.data(), 2*p_-2, nn_);
    ntter_w->ntt_batch(buf.data(), 2*nn_);
    transpose_rect_restrict(dst.data(), buf.data(), nn_, 2*p_-2);
}
void U64Context::iw_intt(vec64& dst, const vec64& src) const
{
    // W-iNTT
    vec64& buf = buf_size_;
    transpose_rect_restrict(buf.data(), src.data(), 2*p_-2, nn_);
    ntter_w->intt_batch(buf.data(), 2*nn_);
    transpose_rect_restrict(dst.data(), buf.data(), nn_, 2*p_-2);
    // I-iNTT
    for(int i=0,j=pnn_;i<pnn_;i++,j++)
    {
        uint64_t Pi = dst[i], Ni = dst[j];
        dst[i] = mod_add(Pi, Ni, q_);
        dst[j] = mm_.mul(
            mod_sub(Pi, Ni, q_),
            I_inv_mont_
        );
    }
    // 除以2
    for(auto& i:dst)i = mm_.mul(i, inv2_mont);
}
void U64Context::xy_ntt(vec64& dst, const vec64& src) const
{
    // x-ntt
    const size_t n = this->n_;
    const size_t nn = this->nn_;
    const size_t pnn = this->pnn_;
    const size_t pn = n*(p_-1);
    // 转置
    for(int i=0;i<2*pnn;i+=nn)transpose_auto(dst.data()+i, src.data()+i, n);
    ntter_p->ntt_batch(dst.data(), pn);
    ntter_n->ntt_batch(dst.data() + pnn, pn);
    // 再转置
    for(int i=0;i<2*pnn;i+=nn)transpose_inplace(dst.data()+i, n);
    ntter_n->ntt_batch(dst.data(), pn);
    ntter_p->ntt_batch(dst.data() + pnn, pn);

}
void U64Context::xy_intt(vec64& dst, const vec64& src) const
{
    // x-ntt
    const size_t n = this->n_;
    const size_t nn = this->nn_;
    const size_t pnn = this->pnn_;
    const size_t pn = n*(p_-1);
    // 转置
    for(int i=0;i<2*pnn;i+=nn)transpose_auto(dst.data()+i, src.data()+i, n);
    ntter_p->intt_batch(dst.data(), pn);
    ntter_n->intt_batch(dst.data() + pnn, pn);
    // 再转置
    for(int i=0;i<2*pnn;i+=nn)transpose_inplace(dst.data()+i, n);
    ntter_n->intt_batch(dst.data(), pn);
    ntter_p->intt_batch(dst.data() + pnn, pn);
}
// 逐位加法
void U64Context::add(vec64& dst, const vec64& src1, const vec64& src2) const
{
    for(int i=0;i<size_;i++)
    {
        dst[i] = mod_add(src1[i], src2[i], q_);
    }
}
// 逐位减法
void U64Context::sub(vec64& dst, const vec64& src1, const vec64& src2) const
{
    for(int i=0;i<size_;i++)
    {
        dst[i] = mod_sub(src1[i], src2[i], q_);
    }
}
void U64Context::sub_unsafe(vec64& dst, const vec64& src1, const vec64& src2) const
{
    for(int i=0;i<size_;i++)
    {
        dst[i] = mod_sub(src1[i] % q_, src2[i] % q_, q_);
    }
}
// 逐位乘法
void U64Context::mul(vec64& dst, const vec64& src1, const vec64& src2) const
{
    for(int i=0;i<size_;i++)
    {
        dst[i] = mod_mul(src1[i], src2[i], q_);
    }
}
// 逐位负
void U64Context::neg(vec64& dst, const vec64& src1) const
{
    for(int i=0;i<size_;i++)
    {
        dst[i] = (q_ - src1[i]) % q_;
    }
}
// 标量乘
void U64Context::mul_scalar(vec64& dst, const vec64& src_vec, uint64_t src_scalar) const
{
    uint64_t encoded_scalar = mm_.encode(src_scalar);
    size_t size = dst.size();
    assert(src_vec.size() == size);
    for(int i=0; i<size; i++)dst[i] = mm_.mul(src_vec[i], encoded_scalar);
}
// 比较
bool U64Context::eq(const vec64& src1, const vec64& src2) const
{
    for(int i=0;i<size_;i++)
    {
        if (src1[i] != src2[i])return false;
    }
    return true;
}