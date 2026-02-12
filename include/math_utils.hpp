#pragma once
#include <cassert>
#include <cstddef>
#include <cstdint>

// 判断一个整数是不是power of 2
bool is_power_of_two(size_t x) noexcept;

size_t Log2(size_t x);

// 一个数在模2^64意义下的乘法逆元
uint64_t modinv64(uint64_t a);