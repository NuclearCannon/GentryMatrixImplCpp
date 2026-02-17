#pragma once
#include <cstdint>
#include <cassert>

// Montgomery乘法器
class MontgomeryMultiplier {
private:
    // 设置私有辅助函数以避免对外依赖


    inline static uint64_t _get_R2(uint64_t M)
    {
        // 其实就是2^128 % M。需要7次平方
        // 使用128位整数以避免乘法溢出
        __uint128_t x = 2;
        x = (x*x) % M;
        x = (x*x) % M;
        x = (x*x) % M;
        x = (x*x) % M;
        x = (x*x) % M;
        x = (x*x) % M;
        x = (x*x) % M;
        return (uint64_t)x;
    }

    inline static uint64_t _modinv64(uint64_t a) {
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


public:
    // M: 模数
    // R: 2^64 % M (在这里隐去)
    // R2: R^2 % M
    // N1: - M^{-1} % 2^{64}
    const uint64_t M, R2, N1;

    MontgomeryMultiplier(uint64_t M):
        M(M),
        R2(_get_R2(M)),
        N1(-_modinv64(M))
    {
        assert(M % 2 == 1);
    }

    // 逻辑上，这相当于return t*R^{-1} mod M
    inline __attribute__((always_inline)) uint64_t montgomery_reduce(__uint128_t t) const {
        uint64_t m = (uint64_t)t * N1;
        uint64_t v = ((t + (__uint128_t)m * M) >> 64);
        return (v<M)?v:v-M;
    }
    
    inline __attribute__((always_inline)) uint64_t mul(uint64_t a, uint64_t b) const {
        return montgomery_reduce((__uint128_t)a * b);
    }

    // 将普通形式整数编码为Montgomery形式
    uint64_t encode(uint64_t x) const {
        return montgomery_reduce((__uint128_t)x * R2);
    }

    // 将Montgomery形式整数还原为普通形式
    uint64_t decode(uint64_t x) const {
        return montgomery_reduce(x);
    }

    
};