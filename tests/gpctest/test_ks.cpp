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