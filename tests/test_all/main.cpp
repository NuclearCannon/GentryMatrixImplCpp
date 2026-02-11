#include "FHE/encrypt64.hpp"
#include "FHE/key_switch_64.hpp"
#include "FHE/circledast.hpp"
#include <chrono>
#include <gperftools/profiler.h>
#include "complecx_matrix.hpp"
#include <iostream>


int main()
{
    // 准备参数
    int n = 256;
    int p = 17;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    u64 qo = 576460752303421441;
    u64 qor = 19;
    std::shared_ptr<U64CtxChain> cc = std::make_shared<U64CtxChain>(n, p, mods, roots);
    double delta = 10000000;
    printf("生成数据\n");
    ComplexMatrixGroup U = ComplexMatrixGroup::random(20, n, p);
    ComplexMatrixGroup V = ComplexMatrixGroup::random(20, n, p);
    ComplexMatrixGroup UV = U.matmul_ABT(V);

    printf("构造明文\n");
    auto u = CRTArray::from_fmpz_vector(U.encode().to_fmpz_vector(delta), cc);
    auto v = CRTArray::from_fmpz_vector(V.encode().to_fmpz_vector(delta), cc);
    // 构造私钥
    printf("构造私钥\n");
    auto sk = CRTArray::sk(cc);
    printf("构造KSK\n");
    auto [ksk1, ksk2] = create_ksks_for_circledast_ct(sk, qo, qor);
    printf("构造密文\n");

    auto [ua, ub] = encrypt64_no_e(u, sk);
    auto [va, vb] = encrypt64_no_e(v, sk);
    printf("执行运算\n");
    auto [wa, wb] = circledast_ct(ua, ub, va, vb, ksk1, ksk2);
    printf("检查结果\n");
    auto w = decrypt64(wa, wb, sk);
    w.mul_scalar_e(n);
    auto W = ComplexMatrixGroup::from_fmpz_vector(
        w.to_fmpz_vector_centered(), delta*delta, n, p
    ).decode();
    double maxr = 0;
    for(int w=0; w<p-1; w++)
    {
        for(int x=0; x<n; x++)
        {
            for(int y=0; y<n; y++)
            {
                // 我们称之为偏差率
                double r = std::abs((W.at(w,x,y) / UV.at(w,x,y))-complex(1));
                if (r > maxr)maxr = r;
            }
        }
    }
    std::cout << "maxr=" <<maxr << std::endl;
    return 0;

}