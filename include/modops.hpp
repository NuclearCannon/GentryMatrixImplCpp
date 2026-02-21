#pragma once
#include <cstdint>


inline __attribute__((always_inline))
uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t mod) {
    // assert(a<mod);
    // assert(b<mod);
    __uint128_t product = (__uint128_t)a * b;
    return (uint64_t)(product % mod);
}


inline __attribute__((always_inline))
uint64_t mod_add(uint64_t a, uint64_t b, uint64_t mod) {
    uint64_t c = a+b;
    return (c<mod)?c:c-mod;
}


inline __attribute__((always_inline))
uint64_t mod_sub(uint64_t a, uint64_t b, uint64_t mod) {
    uint64_t c = a+mod-b;
    return (c<mod)?c:c-mod;
}


// uint64模幂
uint64_t mod_pow(uint64_t base, uint64_t e, uint64_t mod);


// uint64乘法逆元（通过return x^{mod-2}实现）
uint64_t mod_inv(uint64_t x, uint64_t mod);

// 将x取模到[-q/2,q/2]范围内
inline int64_t mod_centered(uint64_t x, uint64_t q)
{
    x %= q;    // 现在属于[0,q)
    if(x*2>q)x-=q;
    return (int64_t)(x);
}