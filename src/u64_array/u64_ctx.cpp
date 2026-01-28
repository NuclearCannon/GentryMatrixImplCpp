#include "u64_array.hpp"


U64Context::U64Context(int n, int p, u64 q, u64 root_q):
    n_(n),
    p_(p),
    nn_(n*n),
    pnn_((p-1)*n*n),
    size_(2*(p-1)*n*n),
    q_(q)
{
    // 计算zeta. zeta是4n阶本原单位根
    assert((q-1)%((u64)(4*n)) == 0);
    u64 zeta = mod_pow(root_q, (q-1)/((u64)(4*n)), q);
    u64 zeta_inv = mod_inv(zeta, q);
    // 计算eta. eta是p阶本原单位根
    assert((q-1)%((u64)p) == 0);
    u64 eta = mod_pow(root_q, (q-1)/((u64)p), q);
    // 生成各ntter对象
    ntter_p = std::make_unique<TwistedNtterXY64>(n,q,zeta);
    ntter_n = std::make_unique<TwistedNtterXY64>(n,q,zeta_inv);
    ntter_w = std::make_unique<TwistedNtterW64>(p,q,eta);
    // 计算I及其逆元
    I_ = mod_pow(zeta, n, q);
    I_inv_ = mod_pow(I_, 3, q);
    assert(mod_mul(I_, I_inv_, q) == 1);

}

U64Context::~U64Context()
{
    // do nothing
}

void U64Context::iw_ntt(vec64& dst, const vec64& src) const
{
    // ntt-I
    // src的低半部分是real，高半部分是imag
    // result的低半部分是P，高半部分是N
    // P = real + imag*I
    // N = real - imag*I
    vec64 real = copy_from(src, 0, pnn_);
    vec64 imag = copy_from(src, pnn_, pnn_);
    vec_scalar_mul(imag, imag, I_, q_);
    for(int i=0;i<pnn_;i++)
    {
        dst[i] = (real[i] + imag[i]) % q_;
        dst[i+pnn_] = (real[i] + q_ - imag[i]) % q_;
    }
    // 现在dst是I-ntted
    vec64 buf_p(p_-1);
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
                ntter_w->ntt(buf_p, buf_p);
                for(int w=0;w<p_-1;w++)dst[base_y + w*nn_] = buf_p[w];
            }
        }
    }

}
void U64Context::iw_intt(vec64& dst, const vec64& src) const
{
    // W-iNTT
    vec64 buf_p(p_-1);
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
                ntter_w->intt(buf_p, buf_p);
                for(int w=0;w<p_-1;w++)dst[base_y + w*nn_] = buf_p[w];
            }
        }
    }
    // I-iNTT
    vec64 real_add_imag_I = copy_from(dst, 0, pnn_);
    vec64 real_sub_imag_I = copy_from(dst, pnn_, pnn_);
    for(int i=0;i<pnn_;i++)
    {
        dst[i] = real_add_imag_I[i] + real_sub_imag_I[i];
        dst[i+pnn_] = mod_mul(
            real_add_imag_I[i] + q_ - real_sub_imag_I[i],
            I_inv_,
            q_
        );
    }
    // 除以2
    vec_scalar_mul(dst, dst, mod_inv(2,q_), q_);
}
void U64Context::xy_ntt(vec64& dst, const vec64& src) const
{
    // x-ntt
    size_t nn = this->nn_;
    vec64 buf_n(n_);

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
                if (i == 0)ntter_p->ntt(buf_n, buf_n);
                else       ntter_n->ntt(buf_n, buf_n);
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
                if (i == 1)ntter_p->ntt(buf_n, buf_n);
                else       ntter_n->ntt(buf_n, buf_n);
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
    vec64 buf_n(n_);

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
                if (i == 0)ntter_p->intt(buf_n, buf_n);
                else       ntter_n->intt(buf_n, buf_n);
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
                if (i == 1)ntter_p->intt(buf_n, buf_n);
                else       ntter_n->intt(buf_n, buf_n);
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
        dst[i] = (src1[i] + src2[i]) % q_;
    }
}
// 逐位减法
void U64Context::sub(vec64& dst, const vec64& src1, const vec64& src2) const
{
    for(int i=0;i<size_;i++)
    {
        dst[i] = (src1[i] + q_ - src2[i]) % q_;
        dst[i] %= q_;
        while(dst[i] < 0) dst[i]+=q_;
        dst[i] %= q_;
    }
}
// 逐位乘法
void U64Context::mul(vec64& dst, const vec64& src1, const vec64& src2) const
{
    vec_mul(dst, src1, src2, q_);
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
    vec_scalar_mul(dst, src_vec, src_scalar, q_);
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