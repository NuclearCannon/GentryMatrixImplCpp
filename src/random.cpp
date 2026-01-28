#include <random>
#include <iostream>
#include <cmath>
#include "random.hpp"
// TODO: 不安全！我们没有使用密码学安全的随机数发生器，而且离散高斯的取法中使用了不安全的浮点数！
// 浮点计算不仅精度有限，还可能引入侧信道攻击风险

namespace ramdom_generators {

// 全局或局部的随机数引擎（推荐使用 thread_local 避免多线程问题）
thread_local static std::mt19937 rng{std::random_device{}()};

double random_real() {
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

int randint(int a, int b) {
    std::uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

u64 randu64(u64 a, u64 b)
{
    std::uniform_int_distribution<u64> dist(a, b);
    return dist(rng);
}

i64 randi64(i64 a, i64 b)
{
    std::uniform_int_distribution<i64> dist(a, b);
    return dist(rng);
}

// 离散高斯分布的似然，未归一化
static double dg_p(int x, double sigma)
{
    // exp(-x^2/(2*sigma^2)
    double x_f = x;
    return std::exp(-(x_f*x_f)/(2*sigma*sigma));
}

// 简单起见，我们使用离散高斯分布实现所有的噪声和私钥
int dg(int bound, double sigma)
{
    while(1)
    {
        int x = randint(-bound, bound);
        double u = random_real();
        if (u <= dg_p(x, sigma)) return x;
    }
}

}; // ramdom_generators