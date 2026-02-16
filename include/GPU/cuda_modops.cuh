#pragma once
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cstdint>

// 用于将一个[0,2M)范围内的数，取模进[0,M)范围
// 无分支，无模运算
static inline __device__ uint64_t _mod(uint64_t x, uint64_t M)
{
    // 先减一个
    int64_t xi = (int64_t)(x-M);
    // 如果x<M, 那么xi<0, 符号位为1
    // 我们将其左移63位以获取符号位
    int64_t mask = xi >> 63;
    // if xi<0: mask=0b1111...
    // else:    mask=0b0000...
    xi += mask & M; // 仅当xi<0时加M，否则+0
    return (uint64_t)xi;
}

static inline __device__ uint64_t _mod_add(uint64_t a, uint64_t b, uint64_t M)
{
    return _mod(a+b, M);
}

static inline __device__ uint64_t _mod_sub(uint64_t a, uint64_t b, uint64_t M)
{
    // 将_mod手动内联进来
    int64_t xi = a-b;   // (-M,M)
    int64_t mask = xi >> 63;
    xi += mask & M;
    return (uint64_t)xi;
}

// 相当于return t*Rinv mod M
static inline __device__ uint64_t _montgomery_reduce_cuda(__uint128_t t, uint64_t M, uint64_t N1) {
    uint64_t m = (uint64_t)t * N1;
    uint64_t v = ((t + (__uint128_t)m * M) >> 64);
    return _mod(v, M);
}

static inline __device__ uint64_t _mul_cuda(uint64_t a, uint64_t b, uint64_t M, uint64_t N1) {
    return _montgomery_reduce_cuda((__uint128_t)a * b, M, N1);
}

struct ModAdd {
    __device__ uint64_t operator()(uint64_t a, uint64_t b, uint64_t M) const {
        return _mod_add(a, b, M);
    }
};

struct ModSub {
    __device__ uint64_t operator()(uint64_t a, uint64_t b, uint64_t M) const {
        return _mod_sub(a, b, M); // 假设你已有 _mod_sub
    }
};