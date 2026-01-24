#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"




void rader(
    int p, int g,
    fmpz_mod_ctx_t q_ctx,
    const std::vector<int>& gpp,
    const std::vector<int>& gnp,
    fmpz_mod_poly_t a, 
    const fmpz_mod_poly_t b, 
    fmpz_mod_poly_t c,
    const fmpz_mod_poly_t xn1, 
    const fmpz_vector& input,
    fmpz_vector& output
)
{
    // 准备poly a
    // a[i] = x[gnp[i]]
    fmpz_mod_poly_set_ui(a, 0, q_ctx);  // 先清零
    for(int i=0;i<p-1;i++)
    {
        fmpz_mod_poly_set_coeff_fmpz(a,i,input[gnp[i]],q_ctx);
    }
    // 和b相乘得到c，并且取模xn1 = X^{p-1}-1
    // fmpz_mod_poly_mulmod(c,a,b,xn1,q_ctx);

    fmpz_mod_poly_mul(c, a, b, q_ctx);

    // 写回
    fmpz_t tmp;
    fmpz_t tmp2;
    fmpz_init(tmp);
    fmpz_init(tmp2);
    for(int i=0;i<p-1;i++)
    {
        // 读取c的X^i系数
        fmpz_mod_poly_get_coeff_fmpz(
            tmp, c, i, q_ctx
        );
        // 读取c的X^(i+p-1)系数，并加回去
        fmpz_mod_poly_get_coeff_fmpz(
            tmp2, c, i+p-1, q_ctx
        );
        fmpz_add(tmp, tmp, tmp2);
        // output[g^i] = x0+c[i]
        fmpz_mod_add(
            output[gpp[i]], tmp, input[0], q_ctx
        );
    }
    // 计算X0 = sum(input)
    fmpz_set_ui(tmp, 0);
    for(int i=0;i<p;i++)
    {
        fmpz_mod_add_fmpz(tmp, tmp, input[i], q_ctx);
    }
    fmpz_set(output[0], tmp);


    fmpz_clear(tmp);
    fmpz_clear(tmp2);

}


TwistedNtterW::TwistedNtterW(int p,int g,const fmpz_t q,const fmpz_t eta):
    p_(p),g_(g),gpp_(p),gnp_(p),buf1(p), buf2(p)
{
    fmpz_init_set(q_, q);
    fmpz_init(p_inv_);
    fmpz_init_set(eta_, eta);
    fmpz_mod_ctx_init(q_ctx_, q_);
    fmpz_mod_poly_init(a_poly_, q_ctx_);
    fmpz_mod_poly_init(b_poly_, q_ctx_);
    fmpz_mod_poly_init(b2_poly_, q_ctx_);
    fmpz_mod_poly_init(c_poly_, q_ctx_);
    fmpz_mod_poly_init(xn1_, q_ctx_);

    // 计算gpp
    gpp_[0] = 1;
    for(int i=1;i<p;i++)gpp_[i] = (gpp_[i-1] * g) % p;
    // 计算gnp。考虑到g^{p-1} = 1, 我们可以认为……
    for(int i=0;i<p;i++)gnp_[i] = gpp_[(p-1)-i];
    // 计算b. b[i] = eta^(gpp[i])
    fmpz_t tmp;
    fmpz_init(tmp);
    for(int i=0;i<p-1;i++)
    {
        // let tmp = eta^gpp[i]
        fmpz_mod_pow_ui(tmp, eta_, gpp_[i], q_ctx_);
        fmpz_mod_poly_set_coeff_fmpz(
            b_poly_, i, tmp, q_ctx_
        );
    }
    // 计算b2
    for(int i=0;i<p-1;i++)
    {
        // let tmp = eta^(-gpp[i])
        fmpz_mod_pow_ui(tmp, eta_, p-gpp_[i], q_ctx_);
        fmpz_mod_poly_set_coeff_fmpz(
            b2_poly_, i, tmp, q_ctx_
        );
    }
    // 计算p_inv_
    fmpz_set_ui(tmp, p_);
    fmpz_invmod(p_inv_, tmp, q_);
    // 设置xn1
    fmpz_mod_poly_set_ui(xn1_, 0, q_ctx_);  // 置为0
    fmpz_mod_poly_set_coeff_si(xn1_, 0, -1, q_ctx_);
    fmpz_mod_poly_set_coeff_si(xn1_, p-1, 1, q_ctx_);


    // finish
    fmpz_clear(tmp);
}
TwistedNtterW::~TwistedNtterW()
{
    fmpz_clear(q_);
    fmpz_clear(eta_);
    fmpz_clear(p_inv_);
    fmpz_mod_ctx_clear(q_ctx_);
    fmpz_mod_poly_clear(a_poly_, q_ctx_);
    fmpz_mod_poly_clear(b_poly_, q_ctx_);
    fmpz_mod_poly_clear(b2_poly_, q_ctx_);
    fmpz_mod_poly_clear(c_poly_, q_ctx_);
    fmpz_mod_poly_clear(xn1_, q_ctx_);
}

void TwistedNtterW::ntt(const fmpz_vector& src, fmpz_vector& dst)
{
    // buf1 = src + [0]
    _fmpz_vec_set(buf1.raw(), src.raw(), p_-1);
    fmpz_set_ui(buf1[p_-1], 0);
    // buf2 = rader(buf1)
    rader(
        p_, g_, q_ctx_, gpp_, gnp_, 
        a_poly_, b_poly_, c_poly_, xn1_, 
        buf1, buf2
    );
    // 写回
    _fmpz_vec_set(dst.raw(), buf2.raw()+1, p_-1);


}
void TwistedNtterW::intt(const fmpz_vector& src, fmpz_vector& dst)
{
    // buf1 = [0]+src
    _fmpz_vec_set(buf1.raw()+1, src.raw(), p_-1);
    fmpz_set_ui(buf1[0], 0);
    // iNTT rader
    rader(
        p_, g_, q_ctx_, gpp_, gnp_, 
        a_poly_, b2_poly_, c_poly_, xn1_, 
        buf1, buf2
    );  // 只要使用b2_poly就好了
    // delta=buf2[p-1]
    for(int i=0;i<p_-1;i++)
    {
        fmpz_mod_sub_fmpz(
            buf1[i],
            buf2[i],
            buf2[p_-1],
            q_ctx_
        );
    }
    // 对buf1除以p，输出到dst
    for(int i=0;i<p_-1;i++)
    {
        fmpz_mod_mul_fmpz(
            dst[i],
            buf1[i],
            p_inv_,
            q_ctx_
        );
    }
    // finish
}