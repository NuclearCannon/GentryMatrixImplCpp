#include "u64_array.hpp"


int test_cuda_add()
{
    printf("测试cuda加法...");
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
    auto u = CRTArray::randint(cc, 1000, 2000);
    auto v = CRTArray::randint(cc, 1000, 2000);
    auto uv_cpu = u.add(v);
    CRTArrayGPU u_gpu(cc), v_gpu(cc), uv_gpu(cc);
    u_gpu.set_from_vv64(u.get_data());
    v_gpu.set_from_vv64(v.get_data());
    uv_gpu.add(u_gpu, v_gpu);

    auto result_gpu = CRTArray(uv_gpu.export_to_vv64(), cc);
    assert(uv_cpu.eq(result_gpu));
    printf("cuda加法测试通过\n");
    return 1;

}


int main()
{
    test_cuda_add();
}