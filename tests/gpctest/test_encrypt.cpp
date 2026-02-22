#include "FHE/encrypt_gp.hpp"
#include "CRT.hpp"

void test_encrypt()
{
    // 准备参数
    int n = 16;
    int p = 5;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};

    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
    };

    GentryPolyCtx ctx(n, p, qrp);
    GentryPoly m = GentryPoly::randint(n, p, mods, -0, +0);
    GentryPoly sk = GentryPoly::sk(n, p, mods);
    // 加密
    auto [cta, ctb] = encrypt_gp(m, sk, ctx);
    // 解密
    auto res = decrypt_gp(cta, ctb, sk, ctx);
    auto error = res;
    GentryPoly::sub(error, error, m);
    int64_t error_abs = error.abs();
    printf("test_encrypt: error_abs=%ld\n", error_abs);

}


void test_encrypt_cuda()
{
    // 准备参数
    int n = 16;
    int p = 5;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};

    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
    };

    GentryPolyCtx ctx(n, p, qrp);
    GentryPoly m = GentryPoly::zeros(n, p, mods, GPDevice::CUDA);
    GentryPoly sk = GentryPoly::sk(n, p, mods).to_cuda();
    // 加密
    auto [cta, ctb] = encrypt_gp(m, sk, ctx);
    // 解密
    auto res = decrypt_gp(cta, ctb, sk, ctx);   // 奇怪，m.hash在这一行发生了更改
    auto error = res;
    GentryPoly::sub(error, error, m);
    int64_t error_abs = error.abs();
    printf("test_encrypt_cuda: error_abs=%ld\n", error_abs);

}