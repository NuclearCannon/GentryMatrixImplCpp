#pragma once
#include <sys/types.h>
#include <vector>

typedef u_int64_t u64;
static_assert(sizeof(u64) == 8);
typedef std::vector<u64> vec64;

// 下面是一些计算工具函数


// uint64模乘。调试期间暂不考虑内联它。
u64 mod_mul(u64 a, u64 b, u64 mod);

// uint64模幂
u64 mod_pow(u64 base, u64 e, u64 mod);

// uinr64乘法逆元（通过return x^{mod-2}实现）
u64 mod_inv(u64 x, u64 mod);
