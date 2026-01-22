#include "FHE/encrypt.hpp"
#include "ziq_array.hpp"
#include "FHE/key_switch.hpp"

int test_ks()
{
    // 准备参数
    int n = 16;
    int p = 5;
    int g = 3;
    fmpz_scalar q("2000000641"), root_q("19");
    ZiqArrayContext ctx(n, p, g, q.raw(), root_q.raw());
    fmpz_scalar q2("1024000339841"), root_q2("3");
    ZiqArrayContext ctx2(n, p, g, q2.raw(), root_q2.raw());

    // 生成明文
    fmpz_vector msg_v(2*(p-1)*n*n);
    for(int i=0;i<2*(p-1)*n*n;i++)fmpz_set_si(msg_v[i], 10000*i);
    ZiqArray msg(msg_v, &ctx);
    // 生成私钥
    ZiqArray sk = ctx.dg();
    ZiqArray sk2 = ctx.dg();
    // 加密
    auto [a,b] = encrypt(msg, sk);
    // 制造KSK
    fmpz_scalar B(16);
    KeySwitchingKey ksk(sk, sk2, &ctx, &ctx2, B.raw(), 10);
    // KS sk->sk2
    auto [a2,b2] = ksk.key_switch_big_2(a, b);
    // 解密
    auto result = decrypt(a2, b2, sk2);
    // 检查结果正确性
    fmpz_vector error = result.sub(msg).data().mod_centered(q.raw());
    int error_max = error.max_abs();
    printf("test_ks: error_max=%d\n", error_max);
    error.print();
    // TODO: 噪声还是太大了，而且，很奇怪的是，它似乎和私钥有关系？
    if (error_max < 300)
    {
        printf("test_ks: pass\n");
        return 1;
    }
    else
    {
        printf("test_ks: 噪声太大\n");
        return 0;
    }

}