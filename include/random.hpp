#pragma once


namespace ramdom_generators {
    double random_real();
    int randint(int a, int b);
    uint64_t randu64(uint64_t a, uint64_t b);
    int64_t randi64(int64_t a, int64_t b);
    // 简单起见，我们使用离散高斯分布实现所有的噪声和私钥
    int dg(int bound, double sigma = 3.19);
}

