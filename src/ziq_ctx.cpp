#include "ziq_array.hpp"
#include "random.hpp"

ZiqArrayContext::ZiqArrayContext(int n, int p, int g, const fmpz_t q, const fmpz_t root_q):
    n_(n), p_(p),
    nn_(n*n), pnn_((p-1)*n*n), size_(2*pnn_),
    ntter_p(nullptr),
    ntter_n(nullptr),
    ntter_w(nullptr),
    mm_ctx(new MatmulContext(n, q)),
    buf_p(new fmpz_vector(p-1)),
    buf_n(new fmpz_vector(n)),
    buf_half(new fmpz_vector(pnn_)),
    buf_size(new fmpz_vector(size_))
{
    // zeta = root_q^((q-1)/4*n)
    // eta  = root_q^((q-1)/p)
    fmpz_scalar zeta, eta, q_sub_1, q_div_4n, q_div_p;
    fmpz_scalar n4(4*n), p_mpz(p);
    fmpz_scalar remainder, quotient, tmp;

    // 赋值q_sub_1 = q-1
    fmpz_sub_si(q_sub_1.raw(), q, 1);
    // 检查q-1是不是能被4n整除
    fmpz_fdiv_qr(quotient.raw(), remainder.raw(), q_sub_1.raw(), n4.raw());
    assert(fmpz_is_zero(remainder.raw()));
    // 此时，商=(q-1)/4n, zeta=root_q^商
    fmpz_powm(zeta.raw(), root_q, quotient.raw(), q);

    // 检查q-1是不是能被p整除
    fmpz_fdiv_qr(quotient.raw(), remainder.raw(), q_sub_1.raw(), p_mpz.raw());
    assert(fmpz_is_zero(remainder.raw()));
    // 此时，商=(q-1)/p, eta=root_q^商
    fmpz_powm(eta.raw(), root_q, quotient.raw(), q);


    fmpz_init_set(q_, q);
    fmpz_init(I_);
    fmpz_init(I_inv_);
    fmpz_mod_ctx_init(q_ctx_, q);
    ntter_p = new TwistedNtterXY(n, q, zeta.raw());

    // let tmp = zeta^{-1}
    fmpz_invmod(tmp.raw(), zeta.raw(), q_);
    ntter_n = new TwistedNtterXY(n, q, tmp.raw());
    ntter_w = new TwistedNtterW(p, g, q, eta.raw());
    // let I = zeta^n
    fmpz_mod_pow_ui(I_, zeta.raw(), n, q_ctx_);
    // let I_inv = I^{-1}
    fmpz_invmod(I_inv_, I_, q_);

}

ZiqArrayContext::~ZiqArrayContext()
{
    fmpz_clear(q_);
    fmpz_clear(I_);
    fmpz_clear(I_inv_);
    fmpz_mod_ctx_clear(q_ctx_);
    delete buf_p, buf_n, buf_half, buf_size, mm_ctx;
}


fmpz_vector ZiqArrayContext::iw_ntt(const fmpz_vector& src) const
{
    fmpz_vector result(size_);
    // ntt-I
    // src的低半部分是real，高半部分是imag
    // result的低半部分是P，高半部分是N
    // P = real + imag*I
    // N = real - imag*I
    fmpz_vector& imag_I = *buf_half;
    const fmpz* real = src.raw();
    const fmpz* imag = src.raw() + pnn_;
    fmpz* P = result.raw();
    fmpz* N = result.raw() + pnn_;
    // let imag_I = imag*I
    _fmpz_mod_vec_scalar_mul_fmpz_mod(
        imag_I.raw(), imag, pnn_, I_, q_ctx_
    );
    _fmpz_mod_vec_add(P, real, imag_I.raw(), pnn_, q_ctx_);
    _fmpz_mod_vec_sub(N, real, imag_I.raw(), pnn_, q_ctx_);
    // 现在result已经是i-ntted了
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
                for(int w=0;w<p_-1;w++)fmpz_set(buf_p->at(w), result[base_y + w*nn_]);
                ntter_w->ntt(*buf_p, *buf_p);
                for(int w=0;w<p_-1;w++)fmpz_set(result[base_y + w*nn_], buf_p->at(w));
            }
        }
    }

    return result;
}
fmpz_vector ZiqArrayContext::iw_intt(const fmpz_vector& src) const
{
    fmpz_vector result(size_);

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
                for(int w=0;w<p_-1;w++)fmpz_set(buf_p->at(w), src[base_y + w*nn_]);
                ntter_w->intt(*buf_p, *buf_p);
                for(int w=0;w<p_-1;w++)fmpz_set(buf_size->at(base_y + w*nn_), buf_p->at(w));
            }
        }
    }
    // 现在buf_size = w-intt(src)
    // 接下来开始i-intt
    fmpz* real = result.raw();
    fmpz* imag = real + pnn_;
    const fmpz* P = buf_size->raw();
    const fmpz* N = P + pnn_;
    // let real = P+N
    _fmpz_mod_vec_add(real, P, N, pnn_, q_ctx_);
    // let buf_half = P-N
    _fmpz_mod_vec_sub(buf_half->raw(), P, N, pnn_, q_ctx_);
    // let imag = buf_half / I
    _fmpz_mod_vec_scalar_mul_fmpz_mod(
        imag, buf_half->raw(), pnn_, I_inv_, q_ctx_
    );
    // result /= 2
    fmpz_t two;
    fmpz_init_set_ui(two, 2);
    _fmpz_mod_vec_scalar_div_fmpz_mod(result.raw(), result.raw(), size_, two, q_ctx_);
    fmpz_clear(two);

    // 完成
    return result;
}

fmpz_vector ZiqArrayContext::xy_ntt(const fmpz_vector& src) const
{
    // x-ntt
    fmpz_vector result(size_);
    size_t nn = this->nn_;

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
                for(int x=0;x<n_;x++)fmpz_set(buf_n->at(x), src[base_y + x*n_]);
                if (i == 0)ntter_p->ntt(*buf_n, *buf_n);
                else       ntter_n->ntt(*buf_n, *buf_n);
                // 写回去
                for(int x=0;x<n_;x++)fmpz_set(buf_size->at(base_y + x*n_), buf_n->at(x));
            }
        }
    }
    // 现在buf_size已经是x-ntt
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
                for(int y=0;y<n_;y++)fmpz_set(buf_n->at(y), buf_size->at(base_x + y));
                if (i == 1)ntter_p->ntt(*buf_n, *buf_n);
                else       ntter_n->ntt(*buf_n, *buf_n);
                // 写回去
                for(int y=0;y<n_;y++)fmpz_set(result[base_x + y], buf_n->at(y));
            }
        }
    }

    return result;
}


fmpz_vector ZiqArrayContext::xy_intt(const fmpz_vector& src) const
{
    // 和上面的区别仅仅在于使用了intt而不是ntt，嗯
    fmpz_vector result(size_);
    size_t nn = this->nn_;

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
                for(int x=0;x<n_;x++)fmpz_set(buf_n->at(x), src[base_y + x*n_]);
                if (i == 0)ntter_p->intt(*buf_n, *buf_n);
                else       ntter_n->intt(*buf_n, *buf_n);
                // 写回去
                for(int x=0;x<n_;x++)fmpz_set(buf_size->at(base_y + x*n_), buf_n->at(x));
            }
        }
    }
    // 现在buf_size已经是x-ntt
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
                for(int y=0;y<n_;y++)fmpz_set(buf_n->at(y), buf_size->at(base_x + y));
                if (i == 1)ntter_p->intt(*buf_n, *buf_n);
                else       ntter_n->intt(*buf_n, *buf_n);
                // 写回去
                for(int y=0;y<n_;y++)fmpz_set(result[base_x + y], buf_n->at(y));
            }
        }
    }

    return result;
}


// 逐位加法
void ZiqArrayContext::add(fmpz_vector& dst, const fmpz_vector& src1, const fmpz_vector& src2) const
{
    int L = dst.len();
    assert(src1.len() == L);
    assert(src2.len() == L);
    _fmpz_mod_vec_add(dst.raw(), src1.raw(), src2.raw(), L, q_ctx_);

}
// 逐位减法
void ZiqArrayContext::sub(fmpz_vector& dst, const fmpz_vector& src1, const fmpz_vector& src2) const
{
    int L = dst.len();
    assert(src1.len() == L);
    assert(src2.len() == L);
    _fmpz_mod_vec_sub(dst.raw(), src1.raw(), src2.raw(), L, q_ctx_);
}
// 逐位乘法
void ZiqArrayContext::mul(fmpz_vector& dst, const fmpz_vector& src1, const fmpz_vector& src2) const
{
    int L = dst.len();
    assert(src1.len() == L);
    assert(src2.len() == L);
    _fmpz_mod_vec_mul(dst.raw(), src1.raw(), src2.raw(), L, q_ctx_);
}
// 逐位负
void ZiqArrayContext::neg(fmpz_vector& dst, const fmpz_vector& src1) const
{
    int L = dst.len();
    assert(src1.len() == L);
    _fmpz_mod_vec_neg(dst.raw(), src1.raw(), L, q_ctx_);
}
// 标量乘
void ZiqArrayContext::mul_scalar(fmpz_vector& dst, const fmpz_vector& src_vec, const fmpz_t src_scalar) const
{
    int L = dst.len();
    assert(src_vec.len() == L);
    _fmpz_mod_vec_scalar_div_fmpz_mod(dst.raw(), src_vec.raw(), L, src_scalar, q_ctx_);
}

bool ZiqArrayContext::eq(const fmpz_vector& src1, const fmpz_vector& src2) const
{
    int L = src1.len();
    assert(src2.len() == L);
    for(int i=0;i<L;i++)
    {
        if(fmpz_cmp(src1[i], src2[i]) != 0)return false;
    }
    return true;
}


ZiqArray ZiqArrayContext::zeros() const
{
    return ZiqArray(fmpz_vector::zeros(size_), this);
}
ZiqArray ZiqArrayContext::uniform() const
{
    fmpz_vector r = fmpz_vector::uniform(size_, q_);
    return ZiqArray(r, this);
}
ZiqArray ZiqArrayContext::randint(int lb, int ub) const
{
    fmpz_vector r = fmpz_vector::zeros(size_);
    for(int i=0;i<size_;i++)fmpz_set_si(r[i], ramdom_generators::randint(lb, ub));
    return ZiqArray(r, this);
}
ZiqArray ZiqArrayContext::dg() const
{
    fmpz_vector r = fmpz_vector::dg(size_);
    // 取个模再说
    // _fmpz_vec_scalar_mod_fmpz(r.raw(), r.raw(), size_, q_);
    return ZiqArray(r, this);
}

ZiqArray ZiqArrayContext::sk() const
{
    fmpz_vector dg = fmpz_vector::dg(2*(p_-1)*n_);
    fmpz_vector r = fmpz_vector::zeros(size_);
    int j = 0;

    for(int i=0;i<2;i++)
    {
        for(int w=0;w<p_-1;w++)
        {
            for(int x=0;x<n_;x++)
            {
                fmpz_set(r[i*pnn_ + w*nn_ + x*n_], dg[j++]);
            }
        }
    }
    
    return ZiqArray(r, this);
}