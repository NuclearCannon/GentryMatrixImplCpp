#pragma once
#include <cassert>
#include <cstddef>
#include <cstdint>

// 判断一个整数是不是power of 2
bool is_power_of_two(size_t x) noexcept;

size_t Log2(size_t x);

// 一个数在模2^64意义下的乘法逆元
uint64_t modinv64(uint64_t a);


void transpose_auto(uint64_t* dst, const uint64_t* src, int n);
void transpose_restrict(uint64_t* __restrict__ dst, const uint64_t* __restrict__ src, int n);
void transpose_inplace(uint64_t* dst, int n);

// 把一个shape=(r, c)的长方形矩阵给转置成(c, r)的
// 行优先编码
void transpose_rect_restrict(uint64_t* __restrict__ dst, const uint64_t* __restrict__ src, size_t r, size_t c);