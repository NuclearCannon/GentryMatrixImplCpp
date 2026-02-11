#include "complecx_matrix.hpp"
#include <cassert>
#include "uint64.hpp"

// exp(i*pi/2n)。zeta是4n阶本原单位根。
complex get_zeta_powers(ssize_t n, ssize_t exp)
{
    return std::polar<double>(1, M_PI_2 * (exp % (4*n)) / n);
}


void naive_dft_XY_complex(
    complex* dst,
    const complex* src,
    ssize_t n,
    bool using_conj_zetas
)
{
    assert(dst != src);
    ssize_t k5 = using_conj_zetas?-1:1;  // k5 = 5^k mod (4n); using_conj_zetas时，通过k5传递“共轭”这一信息
    for(ssize_t k=0; k<n; k++)
    {
        complex Xk(0);
        for (ssize_t j=0; j<n; j++)
        {
            Xk += src[j] * get_zeta_powers(n, k5 * j);
        }
        dst[k] = Xk;
        k5 *= 5;
        k5 %= 4*n;
    }
}

void naive_idft_XY_complex(
    complex* dst,
    const complex* src,
    ssize_t n,
    bool using_conj_zetas
)
{
    assert(dst != src);
    for(ssize_t j=0; j<n; j++)
    {
        complex Xj(0);
        ssize_t k5 = using_conj_zetas?1:-1;
        for (ssize_t k=0; k<n; k++)
        {
            Xj += src[k] * get_zeta_powers(n, k5 * j);
            k5 *= 5;
            k5 %= 4*n;
        }
        dst[j] = Xj / ((complex)(n));
    }
}

// eta是p阶本原单位根，本函数返回eta^exp
complex get_eta_powers(ssize_t p, ssize_t exp)
{
    return std::polar<double>(1, M_PI * 2 * (exp % p) / p);
}

void naive_dft_W_complex(
    complex* dst,
    const complex* src,
    ssize_t p
)
{
    assert(dst != src);
    ssize_t n = p-1;
    assert(is_power_of_two(n) && n>=4);
    // 那么3是一个原根
    constexpr ssize_t gamma = 3;

    ssize_t gamma_k = gamma;
    for(ssize_t k=1; k<p; k++)
    {
        complex Xk = 0;
        for(ssize_t j=0; j<p-1; j++)
        {
            Xk += src[j] * get_eta_powers(p, j*gamma_k);
        }
        dst[k-1] = Xk;
        gamma_k = (gamma_k * gamma) % p;
    }
}

void naive_idft_W_complex(
    complex* dst,
    const complex* src,
    ssize_t p
)
{
    assert(dst != src);
    ssize_t n = p-1;
    assert(is_power_of_two(n) && n>=4);
    // 那么3是一个原根
    constexpr ssize_t gamma = 3;
    
    complex sum = 0;
    for(ssize_t j=0; j<p-1; j++)
    {
        complex Xj = 0;
        ssize_t gamma_k = -gamma;
        for(ssize_t k=1; k<p; k++)
        {
            Xj += src[k-1] * get_eta_powers(p, j*gamma_k);
            gamma_k = (gamma_k * gamma) % p;
        }
        dst[j] = Xj;
        sum += Xj;
    }
    for(int i=0; i<n; i++)
    {
        dst[i] = (dst[i]+sum)/((complex)p);
    }
}

void encode3d(
    complex* dst,
    const complex* src,
    ssize_t n,
    ssize_t p
)
{
    ssize_t buflen = (p-1)<n?n:p-1;
    std::vector<complex> buf1(buflen);
    std::vector<complex> buf2(buflen);

    // W-iDFT
    for(int x=0; x<n; x++)
    {
        for(int y=0; y<n; y++)
        {
            for(int w=0; w<p-1; w++)buf1[w] = src[w*n*n + x*n + y];
            naive_idft_W_complex(buf2.data(), buf1.data(), p);
            for(int w=0; w<p-1; w++)dst[w*n*n + x*n + y] = buf2[w];
        }
    }
    // XY-iDFT
    for(int w=0; w<p-1; w++)
    {
        for(int x=0; x<n; x++)
        {
            for(int y=0; y<n; y++)buf1[y] = dst[w*n*n + x*n + y];
            naive_idft_XY_complex(buf2.data(), buf1.data(), n, true);
            for(int y=0; y<n; y++)dst[w*n*n + x*n + y] = buf2[y];
        }
        for(int y=0; y<n; y++)
        {
            for(int x=0; x<n; x++)buf1[x] = dst[w*n*n + x*n + y];
            naive_idft_XY_complex(buf2.data(), buf1.data(), n, false);
            for(int x=0; x<n; x++)dst[w*n*n + x*n + y] = buf2[x];
        }
    }

}


void decode3d(
    complex* dst,
    const complex* src,
    ssize_t n,
    ssize_t p
)
{
    ssize_t buflen = (p-1)<n?n:p-1;
    std::vector<complex> buf1(buflen);
    std::vector<complex> buf2(buflen);

    
    // XY-DFT
    for(int w=0; w<p-1; w++)
    {
        for(int y=0; y<n; y++)
        {
            for(int x=0; x<n; x++)buf1[x] = src[w*n*n + x*n + y];
            naive_dft_XY_complex(buf2.data(), buf1.data(), n, false);
            for(int x=0; x<n; x++)dst[w*n*n + x*n + y] = buf2[x];
        }
        for(int x=0; x<n; x++)
        {
            for(int y=0; y<n; y++)buf1[y] = dst[w*n*n + x*n + y];
            naive_dft_XY_complex(buf2.data(), buf1.data(), n, true);
            for(int y=0; y<n; y++)dst[w*n*n + x*n + y] = buf2[y];
        }
        
    }
    // W-DFT
    for(int x=0; x<n; x++)
    {
        for(int y=0; y<n; y++)
        {
            for(int w=0; w<p-1; w++)buf1[w] = dst[w*n*n + x*n + y];
            naive_dft_W_complex(buf2.data(), buf1.data(), p);
            for(int w=0; w<p-1; w++)dst[w*n*n + x*n + y] = buf2[w];
        }
    }

}

ComplexMatrixGroup ComplexMatrixGroup::encode() const
{
    ComplexMatrixGroup res(n_, p_);
    encode3d(res.data_.data(), data_.data(), n_, p_);
    return res;

}
ComplexMatrixGroup ComplexMatrixGroup::decode() const
{
    ComplexMatrixGroup res(n_, p_);
    decode3d(res.data_.data(), data_.data(), n_, p_);
    return res;
}