#include "u64_array.hpp"
#include <cstring>


U64Context::U64Context(int n, int p, u64 q, u64 root_q):
    n_(n),
    p_(p),
    nn_(n*n),
    pnn_((p-1)*n*n),
    size_(2*(p-1)*n*n),
    q_(q),
    mm_(q),
    bufp_(p-1),
    bufn_(n)
{
    // 计算zeta. zeta是4n阶本原单位根
    assert((q-1)%((u64)(4*n)) == 0);;
    // 计算eta. eta是p阶本原单位根
    assert((q-1)%((u64)p) == 0);
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
    for(int i=0;i<pnn_;i++)
    {
        u64 real_i = src[i];
        u64 image_I_i = mm_.mul(src[i+pnn_], I_mont_);
        dst[i] = mod_add(real_i, image_I_i, q_);
        dst[i+pnn_] = mod_sub(real_i, image_I_i, q_);
    }
    // 现在dst是I-ntted
    vec64& buf_p = bufp_;
    for(int i=0;i<2;i++)
    {
        size_t base_i = i?pnn_:0;
        // 对data_[base:base+pnn]的区域进行W-NTT
        for(int x=0;x<n_;x++)
        {
            size_t base_x = base_i + n_*x;
            for(int y=0;y<n_;y++)
            {
                size_t base_y = base_x + y;
                for(int w=0;w<p_-1;w++)buf_p[w] = dst[base_y + w*nn_];
                ntter_w->ntt_mont(buf_p, buf_p);
                for(int w=0;w<p_-1;w++)dst[base_y + w*nn_] = buf_p[w];
            }
        }
    }
}
void U64Context::iw_intt(vec64& dst, const vec64& src) const
{
    // W-iNTT
    vec64& buf_p = bufp_;
    for(int i=0;i<2;i++)
    {
        size_t base_i = i?pnn_:0;
        // 对data_[base:base+pnn]的区域进行W-NTT
        for(int x=0;x<n_;x++)
        {
            size_t base_x = base_i + n_*x;
            for(int y=0;y<n_;y++)
            {
                size_t base_y = base_x + y;
                for(int w=0;w<p_-1;w++)buf_p[w] = src[base_y + w*nn_];
                ntter_w->intt_mont(buf_p, buf_p);
                for(int w=0;w<p_-1;w++)dst[base_y + w*nn_] = buf_p[w];
            }
        }
    }
    // I-iNTT
    for(int i=0;i<pnn_;i++)
    {
        u64 Pi = dst[i], Ni = dst[i+pnn_];
        dst[i] = mod_add(Pi, Ni, q_);
        dst[i+pnn_] = mm_.mul(
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
    size_t nn = this->nn_;
    vec64& buf_n = bufn_;

    for(int i=0;i<2;i++)
    {
        size_t base_i = i?pnn_:0;
        for(int w=0;w<p_-1;w++)
        {
            size_t base_w = base_i + w*nn;
            for(int y=0;y<n_;y++)
            {
                size_t base_y = base_w + y;
                // 取出一列
                for(int x=0;x<n_;x++)buf_n[x] = src[base_y + x*n_];
                if (i == 0)ntter_p->ntt_mont(buf_n, buf_n);
                else       ntter_n->ntt_mont(buf_n, buf_n);
                // 写回去
                for(int x=0;x<n_;x++)dst[base_y + x*n_] = buf_n[x];
            }
        }
    }
    // 现在dst已经是x-ntt
    for(int i=0;i<2;i++)
    {
        size_t base_i = i?pnn_:0;
        for(int w=0;w<p_-1;w++)
        {
            size_t base_w = base_i + w*nn;
            for(int x=0;x<n_;x++)
            {
                size_t base_x = base_w + x*n_;
                // 取出一列
                for(int y=0;y<n_;y++)buf_n[y] = dst[base_x + y];
                if (i == 1)ntter_p->ntt_mont(buf_n, buf_n);
                else       ntter_n->ntt_mont(buf_n, buf_n);
                // 写回去
                for(int y=0;y<n_;y++)dst[base_x + y] = buf_n[y];
            }
        }
    }
}
void U64Context::xy_intt(vec64& dst, const vec64& src) const
{
    // x-ntt
    size_t nn = this->nn_;
    vec64& buf_n = bufn_;

    for(int i=0;i<2;i++)
    {
        size_t base_i = i?pnn_:0;
        for(int w=0;w<p_-1;w++)
        {
            size_t base_w = base_i + w*nn;
            for(int y=0;y<n_;y++)
            {
                size_t base_y = base_w + y;
                // 取出一列
                for(int x=0;x<n_;x++)buf_n[x] = src[base_y + x*n_];
                if (i == 0)ntter_p->intt_mont(buf_n, buf_n);
                else       ntter_n->intt_mont(buf_n, buf_n);
                // 写回去
                for(int x=0;x<n_;x++)dst[base_y + x*n_] = buf_n[x];
            }
        }
    }
    for(int i=0;i<2;i++)
    {
        size_t base_i = i?pnn_:0;
        for(int w=0;w<p_-1;w++)
        {
            size_t base_w = base_i + w*nn;
            for(int x=0;x<n_;x++)
            {
                size_t base_x = base_w + x*n_;
                // 取出一列
                for(int y=0;y<n_;y++)buf_n[y] = dst[base_x + y];
                if (i == 1)ntter_p->intt_mont(buf_n, buf_n);
                else       ntter_n->intt_mont(buf_n, buf_n);
                // 写回去
                for(int y=0;y<n_;y++)dst[base_x + y] = buf_n[y];
            }
        }
    }
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
void U64Context::mul_scalar(vec64& dst, const vec64& src_vec, u64 src_scalar) const
{
    u64 encoded_scalar = mm_.encode(src_scalar);
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