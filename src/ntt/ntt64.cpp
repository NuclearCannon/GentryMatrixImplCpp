#include "ntt.hpp"
#include "uint64.hpp"
#include "constant_modulo/mod_int.hpp"

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


// NTT standard（常数模数版）
// M: 模数
// Mr: 模数的一个生成元
// n: 长度（必须为power of 2）
// inverse: true表示ntt逆变换（使用相反zeta, 对结果乘以n.inv()）
template<uint64_t M, uint64_t Mr, size_t n, bool inverse>
void ntt_standard_constant_modulo(u64* dst, const u64* src)
{
    static_assert((n & (n-1)) == 0);
    static_assert((M-1)%n==0);
    assert(&dst != &src);
    const auto& rev = get_bit_reverse_table(n);
    for (size_t i = 0; i < n; ++i) {
        dst[i] = src[rev[i]];
    }
    constexpr ModInt<M> zeta = ModInt<M>(Mr).pow((M-1)/n);
    constexpr ModInt<M> root = inverse?(zeta.inv()):zeta;
    size_t m = 1;
    ModInt<M> exponent, w_m, w, u, v;
    while (m < n) {
        exponent = n / (2 * m);
        w_m = root.pow(exponent.val());
        for (size_t i = 0; i < n; i += 2 * m) {
            w = 1;
            for (size_t j = i; j < i + m; ++j) {
                u = dst[j];
                v = dst[j+m]*w;
                dst[j] = (u + v).val();
                dst[j + m] = (u - v).val();
                w *= w_m;
            }
        }
        m *= 2;
    }
    if constexpr (inverse)
    {
        constexpr ModInt<M> ninv = ModInt<M>(n).inv();
        for(int i=0; i<n; i++)dst[i] = (ModInt<M>(dst[i]) * ninv).val();
    }
}

// 实例化一些
template void ntt_standard_constant_modulo<70368747120641,      6,      256, false>(u64* dst, const u64* src);
template void ntt_standard_constant_modulo<70368747294721,      11,     256, false>(u64* dst, const u64* src);
template void ntt_standard_constant_modulo<70368748426241,      6,      256, false>(u64* dst, const u64* src);
template void ntt_standard_constant_modulo<576460752303421441,  19,     256, false>(u64* dst, const u64* src);

template void ntt_standard_constant_modulo<70368747120641,      6,      256, true >(u64* dst, const u64* src);
template void ntt_standard_constant_modulo<70368747294721,      11,     256, true >(u64* dst, const u64* src);
template void ntt_standard_constant_modulo<70368748426241,      6,      256, true >(u64* dst, const u64* src);
template void ntt_standard_constant_modulo<576460752303421441,  19,     256, true >(u64* dst, const u64* src);


void ntt_standard_64_cm(
    u64* dst, 
    const u64* src, 
    size_t n, 
    const u64 mod,
    bool inverse
)
{
    // 只支持已实例化的组合
    if (n == 256) {
        if (!inverse) {
            if (mod == 70368747120641ULL) {
                ntt_standard_constant_modulo<70368747120641ULL, 6, 256, false>(dst, src);
                return;
            }
            if (mod == 70368747294721ULL) {
                ntt_standard_constant_modulo<70368747294721ULL, 11, 256, false>(dst, src);
                return;
            }
            if (mod == 70368748426241ULL) {
                ntt_standard_constant_modulo<70368748426241ULL, 6, 256, false>(dst, src);
                return;
            }
            if (mod == 576460752303421441ULL) {
                ntt_standard_constant_modulo<576460752303421441ULL, 19, 256, false>(dst, src);
                return;
            }
        } else {
            if (mod == 70368747120641ULL) {
                ntt_standard_constant_modulo<70368747120641ULL, 6, 256, true>(dst, src);
                return;
            }
            if (mod == 70368747294721ULL) {
                ntt_standard_constant_modulo<70368747294721ULL, 11, 256, true>(dst, src);
                return;
            }
            if (mod == 70368748426241ULL) {
                ntt_standard_constant_modulo<70368748426241ULL, 6, 256, true>(dst, src);
                return;
            }
            if (mod == 576460752303421441ULL) {
                ntt_standard_constant_modulo<576460752303421441ULL, 19, 256, true>(dst, src);
                return;
            }
        }
    }
    // 没有对应实例，报错
    throw std::runtime_error("ntt_standard_64_cm: No matching constant-modulo NTT instance for the given parameters.");
}
