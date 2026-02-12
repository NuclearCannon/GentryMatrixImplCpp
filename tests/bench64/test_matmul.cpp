#include "FHE/encrypt64.hpp"
#include "FHE/key_switch_64.hpp"
#include <chrono>
#include <gperftools/profiler.h>

int test_matmul()
{
    // 准备参数
    int n = 256;
    int p = 17;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::shared_ptr<U64CtxChain> cc = std::make_shared<U64CtxChain>(n, p, mods, roots);
    // 构造明文
    auto A = CRTArray::randint(cc, 1000, 2000);
    auto B = CRTArray::randint(cc, 1000, 2000);

    
    auto An = A.iw_ntt();
    auto Bn = B.iw_ntt();
    ProfilerStart("matmul.prof");
    auto Cn = A.circledast(B);
    ProfilerStop();
    // auto C = Cn.iw_intt();
    

    return 1;

}