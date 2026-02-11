#pragma once
#include <vector>
#include <type_traits>
#include <cstdint>
#include <cassert>

typedef std::uint64_t u64;
typedef std::int64_t i64;
static_assert(sizeof(u64) == 8);
typedef std::vector<u64> vec64;
typedef std::vector<vec64> vv64;

// 下面是一些计算工具函数

// TODO: 我们应该在性能测试的情况下避免assert

inline __attribute__((always_inline))
u64 mod_mul(u64 a, u64 b, u64 mod) {
    // assert(a<mod);
    // assert(b<mod);
    __uint128_t product = (__uint128_t)a * b;
    return (u64)(product % mod);
}


inline __attribute__((always_inline))
u64 mod_add(u64 a, u64 b, u64 mod) {
    u64 c = a+b;
    return (c<mod)?c:c-mod;
}


inline __attribute__((always_inline))
u64 mod_sub(u64 a, u64 b, u64 mod) {
    u64 c = a+mod-b;
    return (c<mod)?c:c-mod;
}


// uint64模幂
u64 mod_pow(u64 base, u64 e, u64 mod);


// uinr64乘法逆元（通过return x^{mod-2}实现）
u64 mod_inv(u64 x, u64 mod);




// 判断一个整数是不是power of 2
template<typename T>
inline constexpr bool is_power_of_two(T x) noexcept {
    static_assert(std::is_integral_v<T>, "T must be an integral type");
    return ((x>0) && ((x & (x - 1)) == 0));
}

// 逐位乘（允许别名）
void vec_mul(vec64& dst, const vec64& src1, const vec64& src2, u64 mod);

// 返回x的[0,len)次幂组成的向量
void get_powers(vec64& dst, u64 x, size_t len, u64 mod);

// 返回x的[0,len)次幂组成的向量
vec64 get_powers(u64 x, size_t len, u64 mod);


// 不一定有用，但是调试的时候或许很有用。
std::size_t vec64hash(const vec64& v);
std::size_t vv64hash(const vv64& v);
