#pragma once
#include <cstdint>
#include "modops.hpp"
#include "vec64.hpp"



// Montgomery乘法器
class MontgomeryMultiplier {
private:
    uint64_t M_, R_, Rinv_, R2_, N1_;
public:
    MontgomeryMultiplier(uint64_t M);

    inline uint64_t getM() const {return M_;}
    inline uint64_t getN1() const {return N1_;}
    inline uint64_t getR() const {return R_;}
    inline uint64_t getR2() const {return R2_;}

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