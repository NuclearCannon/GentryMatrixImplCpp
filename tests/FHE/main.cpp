#include "FHE/encrypt.hpp"
#include "ziq_array.hpp"

// 质数: 33532481
// 原根: 3
// 的k : 104789
// zeta: 26242510
// eta: 8938557


int main()
{
    // 准备参数
    int n = 16;
    int p = 5;
    int g = 3;
    fmpz_scalar q("33532481"), zeta("26242510"), eta("8938557");
    ZiqArrayContext ctx(
        n, p, g, q.raw(), zeta.raw(), eta.raw()
    );

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
    // fmpz_vector error = msg.sub(result).data;
    error.print();
    // 理论上应该打印出很多绝对值很小的取值
    return 0;
}