#include "FHE/encrypt.hpp"
#include "ziq_array.hpp"

// 质数: 33532481
// 原根: 3



int test_enc()
{
    // 准备参数
    int n = 16;
    int p = 5;
    int g = 3;
    fmpz_scalar q("33532481"), root_q("3");
    ZiqArrayContext ctx(n, p, g, q.raw(), root_q.raw());

    // 生成明文
    fmpz_vector msg_v(2*(p-1)*n*n);
    for(int i=0;i<2*(p-1)*n*n;i++)fmpz_set_si(msg_v[i], -i);
    ZiqArray msg(msg_v, &ctx);
    // 生成私钥
    ZiqArray sk = ctx.dg();
    // 加密
    auto [a,b] = encrypt(msg, sk);
    // 解密
    auto result = decrypt(a, b, sk);
    // 检查结果正确性
    fmpz_vector error = result.sub(msg).data().mod_centered(q.raw());
    int error_max = error.max_abs();
    printf("test_enc: error_max=%d\n", error_max);
    if (error_max < 20)
    {
        printf("test_enc: pass\n");
        return 1;
    }
    else
    {
        printf("test_enc: 噪声太大\n");
        return 0;
    }
    
}