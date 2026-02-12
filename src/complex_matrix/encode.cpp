#include "complecx_matrix.hpp"
#include <cassert>

#include "ntt.hpp"
#include <cstring>
#include <memory>
#include "math_utils.hpp"

const std::vector<complex>& get_roots(size_t logn, bool conj)
{
    assert(logn < 20);
    int conj_int = conj?1:0;

    static std::unique_ptr<std::vector<complex>> backups[20][2];
    std::unique_ptr<std::vector<complex>>& pos = backups[logn][conj_int];
    if(pos.get())return *pos;
    size_t n = 1<<logn;
    std::vector<complex> roots(logn);
    if (conj) for(int i=0; i<logn; i++)roots[i] = std::polar<double>(1, - 2 * M_PI / double(n>>i));
    else      for(int i=0; i<logn; i++)roots[i] = std::polar<double>(1, + 2 * M_PI / double(n>>i));
    pos = std::make_unique<std::vector<complex>>(std::move(roots));
    return *pos;
}


void dft_standard(complex* dst, const complex* src, size_t n, bool conj)
{
    assert(src != nullptr && dst != nullptr);
    assert(src != dst);
    size_t logn = Log2(n);
    const auto& rev = get_bit_reverse_table_by_logn(logn);
    for (size_t i = 0; i < n; ++i) {
        dst[i] = src[rev[i]];
    }
    size_t m = 1;
    const std::vector<complex>& roots = get_roots(logn, conj);
    size_t log_n_div_2m = logn-1;
    while(m<n)
    {
        complex w_m = roots[log_n_div_2m];
        for(size_t i=0; i<n; i+=2*m)
        {
            complex w = 1;
            for(size_t j=i; j<i+m; j++)
            {
                complex u = dst[j], v=dst[j+m]*w;
                dst[j] = u+v;
                dst[j+m] = u-v;
                w *= w_m;
            }
        }
        m<<=1;
        log_n_div_2m --;
    }
}


void dft_XY_complex(
    complex* dst,
    const complex* src,
    ssize_t n,
    bool using_conj_zetas
)
{
    std::vector<complex> temp(n);
    std::vector<complex> temp2(n);
    // let temp = src * zeta_pos_pows

    
    complex zeta = std::polar<double>(1, M_PI_2 * (using_conj_zetas?-1:1) / n);
    complex zetai = 1;
    for(size_t i=0; i<n; i++)
    {
        temp[i] = src[i] * zetai;
        zetai *= zeta;
    }
    // let temp2 = dst_standard(temp)
    dft_standard(temp2.data(), temp.data(), n, using_conj_zetas);
    // 考虑对应关系
    for(size_t k=0, k5=1; k<n; k++, k5=(k5*5)%(4*n))
    {
        dst[k] = temp2[k5/4];
    }
}

void idft_XY_complex(
    complex* dst,
    const complex* src,
    ssize_t n,
    bool using_conj_zetas
)
{
    std::vector<complex> temp(n);
    std::vector<complex> temp2(n);
    // let temp = src * zeta_pos_pows
    // 考虑对应关系
    for(size_t k=0, k5=1; k<n; k++, k5=(k5*5)%(4*n))
    {
        temp2[k5/4] = src[k];
    }
    // iNTT
    dft_standard(temp.data(), temp2.data(), n, !using_conj_zetas);
    complex izeta = std::polar<double>(1, M_PI_2 * (using_conj_zetas?1:-1) / n);
    complex izetai = 1/((double)n); // 在这除n
    for(size_t i=0; i<n; i++)
    {
        dst[i] = temp[i] * izetai;
        izetai *= izeta;
    }
}


std::vector<complex> get_eta_powers(ssize_t p)
{
    std::vector<complex> result(p);
    int gi = 1;
    for(int i=0; i<p-1; i++, gi=(gi*3)%p)
    {
        result[i] = std::polar<double>(1, M_PI * 2 * gi / p);
    }
    assert(gi==1);
    return result;
}

std::vector<complex> get_ieta_powers(ssize_t p)
{
    std::vector<complex> result(p);
    int gi = 1;
    for(int i=0; i<p-1; i++, gi=(gi*3)%p)
    {
        result[i] = std::polar<double>(1, - M_PI * 2 * gi / p);
    }
    assert(gi==1);
    return result;
}


/*

请不要试图使用Rader来加速W轴的DFT
RaderDFT的数值稳定性太差，而且并不快到哪去！

*/


void naive_dft_W_complex(
    complex* dst,
    const complex* src,
    ssize_t p,
    const std::vector<complex>& etas
)
{
    assert(dst != src);
    ssize_t n = p-1;
    assert(is_power_of_two(n) && n>=4);
    for(int u=0; u<p-1; u++)
    {
        complex Gu = 0;
        for(int v=0; v<p-1; v++)
        {
            Gu += src[v] * etas[(u-v+p-1)%(p-1)];
        }
        dst[u] = Gu;
    }

}

void naive_idft_W_complex(
    complex* dst,
    const complex* src,
    ssize_t p,
    const std::vector<complex>& ietas
)
{
    assert(dst != src);
    ssize_t n = p-1;
    assert(is_power_of_two(n) && n>=4);
    complex T = 0;
    for(int v=0; v<p-1; v++)
    {
        complex gv = 0;
        for(int u=0; u<p-1; u++)
        {
            gv += src[u] * ietas[(u-v+p-1)%(p-1)];
        }
        dst[v] = gv;
        T+=gv;
    }
    for(int i=0; i<p-1; i++){
        dst[i] += T;
        dst[i] /= p;
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
    auto ietas = get_ieta_powers(p);
    // W-iDFT
    for(int x=0; x<n; x++)
    {
        for(int y=0; y<n; y++)
        {
            for(int w=0; w<p-1; w++)buf1[w] = src[w*n*n + x*n + y];
            naive_idft_W_complex(buf2.data(), buf1.data(), p, ietas);
            for(int w=0; w<p-1; w++)dst[w*n*n + x*n + y] = buf2[w];
        }
    }
    // XY-iDFT
    for(int w=0; w<p-1; w++)
    {
        for(int x=0; x<n; x++)
        {
            for(int y=0; y<n; y++)buf1[y] = dst[w*n*n + x*n + y];
            idft_XY_complex(buf2.data(), buf1.data(), n, true);
            for(int y=0; y<n; y++)dst[w*n*n + x*n + y] = buf2[y];
        }
        for(int y=0; y<n; y++)
        {
            for(int x=0; x<n; x++)buf1[x] = dst[w*n*n + x*n + y];
            idft_XY_complex(buf2.data(), buf1.data(), n, false);
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
    auto etas = get_eta_powers(p);
    
    // XY-DFT
    for(int w=0; w<p-1; w++)
    {
        for(int y=0; y<n; y++)
        {
            for(int x=0; x<n; x++)buf1[x] = src[w*n*n + x*n + y];
            dft_XY_complex(buf2.data(), buf1.data(), n, false);
            for(int x=0; x<n; x++)dst[w*n*n + x*n + y] = buf2[x];
        }
        for(int x=0; x<n; x++)
        {
            for(int y=0; y<n; y++)buf1[y] = dst[w*n*n + x*n + y];
            dft_XY_complex(buf2.data(), buf1.data(), n, true);
            for(int y=0; y<n; y++)dst[w*n*n + x*n + y] = buf2[y];
        }
        
    }
    // W-DFT
    for(int x=0; x<n; x++)
    {
        for(int y=0; y<n; y++)
        {
            for(int w=0; w<p-1; w++)buf1[w] = dst[w*n*n + x*n + y];
            naive_dft_W_complex(buf2.data(), buf1.data(), p, etas);
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