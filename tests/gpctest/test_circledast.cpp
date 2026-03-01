#include "FHE/encrypt_gp.hpp"
#include "FHE/key_switch_gp.hpp"
#include "FHE/circledast.hpp"
#include "CRT.hpp"

void test_circledast()
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
    GentryPoly u = GentryPoly::randint(n, p, mods, -1000, 1000);
    GentryPoly v = GentryPoly::randint(n, p, mods, -1000, 1000);
    GentryPoly sk = GentryPoly::sk(n, p, mods);

    auto [ua, ub] = encrypt_gp(u, sk, ctx);
    auto [va, vb] = encrypt_gp(v, sk, ctx);

    // 生成含有qo的sk1, sk2
    // 生成KSK
    auto ksk_pair = create_ksks_for_circledast_ct(sk, qo, ctx);
    auto u2 = decrypt_gp(ua, ub, sk, ctx);
    auto v2 = decrypt_gp(va, vb, sk, ctx);

    auto [ra, rb] = circledast_ct(ua, ub, va, vb, ksk_pair.first, ksk_pair.second, ctx);
    GentryPoly r = decrypt_gp(ra, rb, sk, ctx);
    // 直接计算
    GentryPoly w = GentryPoly::zeros_like(u);
    u2.iw_ntt(ctx);
    v2.iw_ntt(ctx);
    GentryPoly::circledast(w, u2, v2);  // 使用u2 @ v2作为对照组，因为我们忽略加密本身的噪声
    w.iw_intt(ctx);
    GentryPoly::sub(r,r,w);
    printf("test_circledast: error_abs=%ld\n", r.abs());

}

