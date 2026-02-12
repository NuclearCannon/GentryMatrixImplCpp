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
class MontgomeryMultiplier {
private:
    uint64_t M_, R_, Rinv_, R2_, N1_;
public:
    MontgomeryMultiplier(uint64_t M);

    inline uint64_t getM() const {return M_;}

    // 相当于return t*Rinv mod M
    inline __attribute__((always_inline)) uint64_t montgomery_reduce(__uint128_t t) const {
        uint64_t m = (uint64_t)t * N1_;
        uint64_t v = ((t + (__uint128_t)m * M_) >> 64);
        return (v<M_)?v:v-M_;
    }
    
    inline __attribute__((always_inline)) uint64_t mul(uint64_t a, uint64_t b) const {
        return montgomery_reduce((__uint128_t)a * b);
    }

    // 将普通形式整数编码为Montgomery形式
    uint64_t encode(uint64_t x) const {
        return montgomery_reduce((__uint128_t)x * R2_);
    }

    // 将Montgomery形式整数还原为普通形式
    uint64_t decode(uint64_t x) const {
        return montgomery_reduce(x);
    }

    vec64 batch_encode(const vec64& src) const;
    void batch_encode_to(vec64& dst, const vec64& src) const;
    void batch_decode_to(vec64& dst, const vec64& src) const;
    void batch_encode_inplace(vec64& v) const;
    void batch_decode_inplace(vec64& v) const;

    void vec_mul_mont(vec64& dst, const vec64& src1, const vec64& src2) const;
    
};