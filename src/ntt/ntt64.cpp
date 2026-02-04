#include "ntt.hpp"
#include "uint64.hpp"
#include "montgomery.hpp"


void ntt_standard_64(
    const u64* a, 
    u64* dst, 
    u64 root,
    size_t n, 
    const u64 mod
)
{
    std::vector<u64> roots(n);
    roots[0] = 1;
    for(int i=1;i<n;i++)roots[i] = mod_mul(roots[i-1], root, mod);
    ntt_standard_64_with_roots(a,dst,roots.data(), n, mod);
}


void ntt_standard_64_with_roots(
    const u64* a, 
    u64* dst, 
    const u64* roots,   // 需要提供root的至少[0,n/2)次方 
    size_t n, 
    const u64 mod
)
{
    assert(a != nullptr && dst != nullptr);
    assert((n & (n - 1)) == 0 && n >= 1);
    assert(a != dst);   // 暂不支持……
    const auto& rev = get_bit_reverse_table(n);
    // dst[i] = a[rev[i]] mod ctx
    for (size_t i = 0; i < n; ++i) {
        dst[i] = a[rev[i]] % mod;
    }
    size_t m = 1;
    int t = log2(n>>1);
    while (t>=0) {
        for (size_t i = 0; i < n; i += 2 * m) {
            for (size_t k=0; k < m; ++k) {
                size_t j = i+k;
                // printf("k=%zu, t=%d, k<<t=%zu\n", k, t, k<<t);
                u64 w = roots[k<<t];
                u64 u = dst[j];
                u64 v = mod_mul(dst[j+m], w, mod);
                dst[j] = mod_add(u, v, mod);
                dst[j+m] = mod_sub(u, v, mod);         
            }
        }
        m <<= 1;    // m *= 2
        t--;
    }
}

template<uint64_t M>
vec64 get_powers_mont(u64 x, int len)
{
    vec64 res = get_powers(x, len, M);
    for(int i=0;i<len;i++)res[i] = MontgomeryMultiplier<M>::encode(res[i]);
    return res;
}
// NTT standard（常数模数版）
// M: 模数
// Mr: 模数的一个生成元
// n: 长度（必须为power of 2）
// inverse: true表示ntt逆变换（使用相反zeta, 对结果乘以n.inv()）
template<uint64_t M, uint64_t Mr, size_t n, bool inverse>
void ntt_standard_constant_modulo(u64* dst, const u64* src, bool mont_in, bool mont_out)
{
    static_assert((n & (n-1)) == 0);
    static_assert((M-1)%n==0);
    assert(&dst != &src);
    constexpr int logn = log2(n);
    static const auto& rev = get_bit_reverse_table_by_logn(logn);
    if (mont_in)
    {
        for (size_t i = 0; i < n; ++i) {
            dst[i] = src[rev[i]];
        }
    }
    else
    {
        for (size_t i = 0; i < n; ++i) {
            dst[i] = MontgomeryMultiplier<M>::encode(src[rev[i]]);
        }
    }
    
    constexpr u64 zeta = mod_pow_tmpl<M>(Mr, (M-1)/n);
    constexpr u64 root = (inverse?(mod_inv_tmpl<M>(zeta)):zeta);
    static const vec64 roots = get_powers_mont<M>(root, n);
    size_t m = 1;
    int t = logn-1;
    while (t>=0) {
        for (size_t i = 0; i < n; i += 2 * m) {
            for (size_t k=0; k < m; ++k) {
                size_t j = i+k;
                u64 w = roots[k<<t];
                u64 u = dst[j];
                u64 v = MontgomeryMultiplier<M>::mul(dst[j+m], w);
                dst[j] = mod_add_tmpl<M>(u, v);
                dst[j+m] = mod_sub_tmpl<M>(u, v);         
            }
        }
        m <<= 1;    // m *= 2
        t--;
    }
    if constexpr (inverse)
    {
        constexpr u64 ninv = MontgomeryMultiplier<M>::encode(mod_inv_tmpl<M>(n));
        for(int i=0; i<n; i++)dst[i] = MontgomeryMultiplier<M>::mul(dst[i], ninv);
    }
    if(!mont_out)
    {
        for(int i=0; i<n; i++)dst[i] = MontgomeryMultiplier<M>::decode(dst[i]);
    }
}

// 实例化一些
template void ntt_standard_constant_modulo<70368747120641,      6,      256, true >(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368747294721,      11,     256, true >(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368748426241,      6,      256, true >(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<576460752303421441,  19,     256, true >(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368747120641,      6,      16, true >(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368747294721,      11,     16, true >(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368748426241,      6,      16, true >(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<576460752303421441,  19,     16, true >(u64* dst, const u64* src, bool, bool);

template void ntt_standard_constant_modulo<70368747120641,      6,      256, false>(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368747294721,      11,     256, false>(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368748426241,      6,      256, false>(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<576460752303421441,  19,     256, false>(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368747120641,      6,      16, false>(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368747294721,      11,     16, false>(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<70368748426241,      6,      16, false>(u64* dst, const u64* src, bool, bool);
template void ntt_standard_constant_modulo<576460752303421441,  19,     16, false>(u64* dst, const u64* src, bool, bool);


void ntt_standard_64_cm(
    u64* dst, 
    const u64* src, 
    size_t n, 
    const u64 mod,
    bool inverse,
    bool mont_in,
    bool mont_out
)
{
    // 只支持已实例化的组合
    if (n == 256) {
        if (!inverse) {
            if (mod == 70368747120641ULL) {
                ntt_standard_constant_modulo<70368747120641ULL, 6, 256, false>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 70368747294721ULL) {
                ntt_standard_constant_modulo<70368747294721ULL, 11, 256, false>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 70368748426241ULL) {
                ntt_standard_constant_modulo<70368748426241ULL, 6, 256, false>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 576460752303421441ULL) {
                ntt_standard_constant_modulo<576460752303421441ULL, 19, 256, false>(dst, src, mont_in, mont_out);
                return;
            }
        } else {
            if (mod == 70368747120641ULL) {
                ntt_standard_constant_modulo<70368747120641ULL, 6, 256, true>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 70368747294721ULL) {
                ntt_standard_constant_modulo<70368747294721ULL, 11, 256, true>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 70368748426241ULL) {
                ntt_standard_constant_modulo<70368748426241ULL, 6, 256, true>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 576460752303421441ULL) {
                ntt_standard_constant_modulo<576460752303421441ULL, 19, 256, true>(dst, src, mont_in, mont_out);
                return;
            }
        }
    }
    if (n == 16) {
        if (!inverse) {
            if (mod == 70368747120641ULL) {
                ntt_standard_constant_modulo<70368747120641ULL, 6, 16, false>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 70368747294721ULL) {
                ntt_standard_constant_modulo<70368747294721ULL, 11, 16, false>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 70368748426241ULL) {
                ntt_standard_constant_modulo<70368748426241ULL, 6, 16, false>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 576460752303421441ULL) {
                ntt_standard_constant_modulo<576460752303421441ULL, 19, 16, false>(dst, src, mont_in, mont_out);
                return;
            }
        } else {
            if (mod == 70368747120641ULL) {
                ntt_standard_constant_modulo<70368747120641ULL, 6, 16, true>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 70368747294721ULL) {
                ntt_standard_constant_modulo<70368747294721ULL, 11, 16, true>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 70368748426241ULL) {
                ntt_standard_constant_modulo<70368748426241ULL, 6, 16, true>(dst, src, mont_in, mont_out);
                return;
            }
            if (mod == 576460752303421441ULL) {
                ntt_standard_constant_modulo<576460752303421441ULL, 19, 16, true>(dst, src, mont_in, mont_out);
                return;
            }
        }
    }
    // 没有对应实例，报错
    throw std::runtime_error(
        "ntt_standard_64_cm: No matching constant-modulo NTT instance for the given parameters. "
        "n=" + std::to_string(n) +
        ", mod=" + std::to_string(mod) +
        ", inverse=" + std::to_string(inverse)
    );
}
