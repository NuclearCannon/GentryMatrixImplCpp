#include "FHE/encrypt.hpp"
#include "ziq_array.hpp"
#include "FHE/key_switch.hpp"
#include "FHE/circledast_ct.hpp"
#include <chrono>

ZiqArray circledast_pt(const ZiqArray& A, const ZiqArray& B)
{
    return A.iw_ntt().circledast(B.iw_ntt()).iw_intt();
}


// q  696898287454081973172991196020261322423297
// gq 3
// Q  803469022129495137770981046170581301261101496891396432432129
// gQ 3
// gp 3

int bench_circledast_ct()
{
    // 准备参数
    int n = 16;
    int p = 17;
    int g = 3;
    fmpz_scalar q("696898287454081973172991196020261322423297"), root_q("3");
    ZiqArrayContext ctx(n, p, g, q.raw(), root_q.raw());
    fmpz_scalar q2("803469022129495137770981046170581301261101496891396432432129"), root_q2("3");
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
    fmpz_scalar B(1<<30);int L=5;

    // 统计 create_ksks_for_circledast 耗时
    auto t1 = std::chrono::high_resolution_clock::now();
    auto [ksk1, ksk2] = create_ksks_for_circledast(sk, &ctx, &ctx2, B.raw(), L);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration_ksk = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    printf("create_ksks_for_circledast: %ld us\n", duration_ksk);

    // 统计 circledast_ct 耗时
    auto t3 = std::chrono::high_resolution_clock::now();
    auto [a3,b3] = circledast_ct(a1,b1,a2,b2,*ksk1,*ksk2);
    auto t4 = std::chrono::high_resolution_clock::now();
    auto duration_circledast = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
    printf("circledast_ct: %ld us\n", duration_circledast);
    ZiqArray msg3 = decrypt(a3,b3,sk);

    // 检查结果正确性
    fmpz_vector error = msg3.sub(msg3_real).data().mod_centered(q.raw());
    int error_max = error.max_abs();
    printf("circledast_ct: error_max=%d\n", error_max);
    return 0;

}