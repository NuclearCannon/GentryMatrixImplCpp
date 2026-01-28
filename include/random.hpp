#pragma once
#include "uint64.hpp"

namespace ramdom_generators {
    double random_real();
    int randint(int a, int b);
    u64 randu64(u64 a, u64 b);
    i64 randi64(i64 a, i64 b);
    // 简单起见，我们使用离散高斯分布实现所有的噪声和私钥
    int dg(int bound, double sigma = 3.19);
}

