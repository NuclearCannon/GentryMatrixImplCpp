#pragma once
#include <cstdint>
#include "uint64.hpp"

// 一个数在模2^64意义下的乘法逆元
constexpr uint64_t modinv64(uint64_t a) {
    // 要求 a 是奇数，否则无逆元
    if ((a & 1) == 0) {
        return 0; 
    }

    uint64_t x = 1; // a ≡ 1 (mod 2), so inverse is 1 mod 2
    // 牛顿迭代：每次加倍精度
    x = x * (2 - a * x); // 2 bits
    x = x * (2 - a * x); // 4 bits
    x = x * (2 - a * x); // 8 bits
    x = x * (2 - a * x); // 16 bits
    x = x * (2 - a * x); // 32 bits
    x = x * (2 - a * x); // 64 bits

    return x;
}


// Montgomery乘法器
template<uint64_t M>
struct MontgomeryMultiplier {
    static_assert(M % 2 == 1, "M must be odd");

    // R = 2^64 mod M，需要预计算
    static constexpr uint64_t R = mod_pow_tmpl<M>(2, 64); 

    // Rinv = pow(R, -1, M) = pow(R, M-2, M)
    static constexpr uint64_t Rinv = mod_pow_tmpl<M>(R, M-2); 
    static_assert(mod_mul_tmpl<M>(Rinv, R) == 1);

    // R2 = R^2 % M
    static constexpr uint64_t R2 = mod_mul_tmpl<M>(R, R);
    // -M^(-1) mod 2^64
    static constexpr uint64_t Minv64 = modinv64(M);
    static_assert(Minv64 * M == 1);

    static constexpr uint64_t N_1 = -Minv64;

    
    // 相当于return t*Rinv mod M
    static constexpr inline __attribute__((always_inline)) uint64_t montgomery_reduce(__uint128_t t) {
        uint64_t m = (uint64_t)t * N_1;
        uint64_t v = ((t + (__uint128_t)m * M) >> 64);
        if(v>=M)v-=M;
        return v;
    }
    
    static constexpr inline __attribute__((always_inline)) uint64_t mul(uint64_t a, uint64_t b) {
        return montgomery_reduce((__uint128_t)a * b);
    }

    // 将普通形式整数编码为Montgomery形式
    static constexpr uint64_t encode(uint64_t x) {
        return montgomery_reduce((__uint128_t)x * R2);
    }

    // 将Montgomery形式整数还原为普通形式
    static constexpr uint64_t decode(uint64_t x) {
        return montgomery_reduce(x);
    }
    
};