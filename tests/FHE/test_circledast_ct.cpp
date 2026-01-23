#include "FHE/encrypt.hpp"
#include "ziq_array.hpp"
#include "FHE/key_switch.hpp"
#include "FHE/circledast_ct.hpp"

ZiqArray circledast_pt(const ZiqArray& A, const ZiqArray& B)
{
    return A.iw_ntt().circledast(B.iw_ntt()).iw_intt();
}


int test_circledast_ct()
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

    ZiqArray msg1 = ctx.randint(-50000, +50000);
    ZiqArray msg2 = ctx.randint(-50000, +50000);
    
    ZiqArray msg3_real = circledast_pt(msg1, msg2);

    // 生成私钥
    ZiqArray sk = ctx.sk();
    // ZiqArray sk = ctx.zeros();
    // 加密
    auto [a1,b1] = encrypt_no_e(msg1, sk);
    auto [a2,b2] = encrypt_no_e(msg2, sk);

    // 制造KSK
    fmpz_scalar B(16);int L=10;

    auto [ksk1, ksk2] = create_ksks_for_circledast(sk, &ctx, &ctx2, B.raw(), L);
    auto [a3,b3] = circledast_ct(a1,b1,a2,b2,*ksk1,*ksk2);
    ZiqArray msg3 = decrypt(a3,b3,sk);

    // 检查结果正确性
    fmpz_vector error = msg3.sub(msg3_real).data().mod_centered(q.raw());
    int error_max = error.max_abs();
    printf("circledast_ct: error_max=%d\n", error_max);
    error.print();
    // TODO: 噪声还是太大了，而且，很奇怪的是，它似乎和私钥有关系？
    if (error_max < 100)
    {
        printf("circledast_ct: pass\n");
        return 1;
    }
    else
    {
        printf("circledast_ct: 噪声太大\n");
        return 0;
    }

}