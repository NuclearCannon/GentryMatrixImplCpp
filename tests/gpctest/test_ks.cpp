#include "FHE/encrypt_gp.hpp"
#include "FHE/key_switch_gp.hpp"
#include "CRT.hpp"

void test_ks()
{
    // 准备参数
    int n = 16;
    int p = 5;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
        {qo, qor}
    };

    GentryPolyCtx ctx(n, p, qrp);
    GentryPoly m = GentryPoly::randint(n, p, mods, -1000, +1000);
    GentryPoly sk1 = GentryPoly::sk(n, p, mods);
    GentryPoly sk2 = GentryPoly::sk(n, p, mods);
    // 加密（使用sk1）
    auto [cta, ctb] = encrypt_gp(m, sk1, ctx);
    // 生成含有qo的sk1, sk2
    GentryPoly sk1qo = sk1, sk2qo = sk2;
    sk1qo.moduli_extend_mult(qo);
    sk2qo.moduli_extend_unsafe(qo);
    // 生成KSK
    KeySwitchKeyGP ksk = KeySwitchKeyGP::ksk_gen(sk1qo, sk2qo, qo, ctx);
    auto [a2, b2]  = ksk.key_switch_big_2(cta, ctb, ctx);
    // 解密
    auto res = decrypt_gp(a2, b2, sk2, ctx);
    auto error = res;
    GentryPoly::sub(error, error, m);
    int64_t error_abs = error.abs();
    printf("test_ks: error_abs=%ld\n", error_abs);

}

void test_ks_cuda()
{
    // 准备参数
    int n = 16;
    int p = 5;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
        {qo, qor}
    };

    GentryPolyCtx ctx(n, p, qrp);
    GentryPoly m = GentryPoly::randint(n, p, mods, -1000, +1000);
    GentryPoly sk1 = GentryPoly::sk(n, p, mods);
    GentryPoly sk2 = GentryPoly::sk(n, p, mods);
    // 加密（使用sk1）
    auto [cta, ctb] = encrypt_gp(m, sk1, ctx);
    auto cta_cuda = cta.to_cuda();
    auto ctb_cuda = ctb.to_cuda();

    // 生成含有qo的sk1, sk2
    GentryPoly sk1qo = sk1, sk2qo = sk2;
    sk1qo.moduli_extend_mult(qo);
    sk2qo.moduli_extend_unsafe(qo);
    // 生成KSK
    KeySwitchKeyGP ksk = KeySwitchKeyGP::ksk_gen(sk1qo, sk2qo, qo, ctx);

    auto [a2, b2]  = ksk.key_switch_big_2(cta, ctb, ctx);
    auto [a3, b3]  = ksk.key_switch_big_2(cta_cuda, ctb_cuda, ctx);
    assert(a2.eq(a3));
    assert(b2.eq(b3));
    // 解密
    auto res_cpu = decrypt_gp(a2, b2, sk2, ctx);
    auto res_gpu = decrypt_gp(a3, b3, sk2.to_cuda(), ctx);
    assert(res_cpu.eq(res_gpu));
    auto error = res_cpu;
    GentryPoly::sub(error, error, m);
    int64_t error_abs = error.abs();
    printf("test_ks_cuda: error_abs=%ld\n", error_abs); // 17

    auto error2 = res_gpu;
    GentryPoly::sub(error2, error2, m.to_cuda());
    int64_t error_abs2 = error2.abs();
    printf("test_ks_cuda: error_abs2=%ld\n", error_abs2);   // -1（表示结果超大）

}