#include "FHE/encrypt.hpp"
#include "ziq_array.hpp"
#include "FHE/key_switch.hpp"

int bench_ks()
{
    // 准备参数
    int n = 256;
    int p = 17;
    int g = 3;
    // fmpz_scalar q("696898287454081973172991196020261322423297"), root_q("3");
    fmpz_scalar q("549774251009"), root_q("3");
    ZiqArrayContext ctx(n, p, g, q.raw(), root_q.raw());
    // fmpz_scalar q2("803469022129495137770981046170581301261101496891396432432129"), root_q2("3");
    fmpz_scalar q2("633825300114114700748416654337"), root_q2("3");
    ZiqArrayContext ctx2(n, p, g, q2.raw(), root_q2.raw());
    // 生成私钥
    ZiqArray sk = ctx.dg();
    ZiqArray sk2 = ctx.sk();
    // 加密
    // 制造KSK
    fmpz_scalar B(1<<30);
    KeySwitchingKey ksk(sk, sk2, &ctx, &ctx2, B.raw(), 2);
    // KS sk->sk2
    ZiqArray a = ctx.uniform();
    auto [a2,b2] = ksk.key_switch_big_1(a);
    return 0;
}