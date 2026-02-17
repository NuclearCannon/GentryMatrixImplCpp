#include "u64_array.hpp"


int test_cuda_add_sub()
{
    printf("测试cuda加减法...");
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
    auto u_add_v = u.add(v);
    auto u_sub_v = u.sub(v);
    CRTArrayGPU u_gpu(cc), v_gpu(cc), w_gpu(cc);
    u_gpu.set_from_vv64(u.get_data());
    v_gpu.set_from_vv64(v.get_data());
    w_gpu.add(u_gpu, v_gpu);
    auto result_gpu = CRTArray(w_gpu.export_to_vv64(), cc);
    assert(u_add_v.eq(result_gpu));
    w_gpu.sub(u_gpu, v_gpu);
    auto result_gpu2 = CRTArray(w_gpu.export_to_vv64(), cc);
    assert(u_sub_v.eq(result_gpu2));

    auto neg_u = u.neg();
    u_gpu.neg_inplace();
    auto result_gpu3 = CRTArray(u_gpu.export_to_vv64(), cc);
    assert(neg_u.eq(result_gpu3));

    printf("cuda加减法测试通过\n");
    return 1;

}

int test_cuda_mul_mont()
{
    printf("测试cuda蒙哥马利乘法...");
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
    auto u_add_v = u.mul_mont(v);
    CRTArrayGPU u_gpu(cc), v_gpu(cc), w_gpu(cc);
    u_gpu.set_from_vv64(u.get_data());
    v_gpu.set_from_vv64(v.get_data());
    w_gpu.mul_mont(u_gpu, v_gpu);
    auto result_gpu = CRTArray(w_gpu.export_to_vv64(), cc);
    assert(u_add_v.eq(result_gpu));
    printf("cuda蒙哥马利乘法测试通过\n");
    return 1;

}


int test_cuda_scalar_mul()
{
    printf("测试cuda标量乘法...");
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
    uint64_t v = 1234;
    auto uv = u.mul_scalar(v);
    CRTArrayGPU u_gpu(cc);
    u_gpu.set_from_vv64(u.get_data());
    u_gpu.mul_scalar_inplace(v);
    assert(uv.eq(CRTArray(u_gpu.export_to_vv64(), cc)));
    printf("cuda标量乘法测试通过\n");
    return 1;
}


int test_cuda_mont_encode_decode()
{
    printf("测试cuda蒙哥马利编解码...");
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
    auto u_encoded = u.mont_encode();
    auto u_decoded = u.mont_decode();
    CRTArrayGPU u_gpu(cc);
    u_gpu.set_from_vv64(u.get_data());
    u_gpu.mont_encode_inplace();
    assert(u_encoded.eq(CRTArray(u_gpu.export_to_vv64(), cc)));
    u_gpu.set_from_vv64(u.get_data());
    u_gpu.mont_decode_inplace();
    assert(u_decoded.eq(CRTArray(u_gpu.export_to_vv64(), cc)));
    printf("cuda蒙哥马利编解码测试通过\n");
    return 1;
}



int test_cuda_iwntt()
{
    printf("测试cuda iw ntt...");
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
    auto u_ntted = u.iw_ntt();
    CRTArrayGPU u_gpu(cc);
    u_gpu.set_from_vv64(u.get_data());
    u_gpu.iw_ntt_inplace();
    assert(u_ntted.eq(CRTArray(u_gpu.export_to_vv64(), cc)));
    printf("cuda iw ntt测试通过\n");
    return 1;
}

int test_cuda_iwintt()
{
    printf("测试cuda iw intt...");
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
    auto u_ntted = u.iw_intt();
    CRTArrayGPU u_gpu(cc);
    u_gpu.set_from_vv64(u.get_data());
    u_gpu.iw_intt_inplace();
    assert(u_ntted.eq(CRTArray(u_gpu.export_to_vv64(), cc)));
    printf("cuda iw intt测试通过\n");
    return 1;
}


int main()
{
    test_cuda_add_sub();
    test_cuda_mul_mont();
    test_cuda_scalar_mul();
    test_cuda_mont_encode_decode();
    test_cuda_iwntt();
    test_cuda_iwintt();
    return 0;
}