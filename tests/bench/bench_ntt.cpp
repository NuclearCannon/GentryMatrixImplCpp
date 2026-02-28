#include "GPU/cuda_ntt.hpp"
#include <gperftools/profiler.h>
#include <vector>
#include "modops.hpp"
#include <iostream>


void bench_ntt_cuda()
{
    // 准备参数
    int logn = 8;
    int n = 1<<logn;    // 单次NTT规模
    int batch_size = 1024;   // 批次大小
    uint64_t q = 70368747120641;// 模数
    uint64_t qr = 6;    // q的原根
    uint64_t omega = mod_pow(qr, (q-1)/n, q);
    uint64_t iomega = mod_inv(omega, q);
    MontgomeryMultiplier mm(q);

    std::cout << "NTT(CUDA)性能测试" << std::endl;
    std::cout << "n=" << n << "\tbatch=" << batch_size << "\tq=" << q << std::endl;

    CudaBuffer buf1(n*batch_size*sizeof(uint64_t));
    CudaBuffer omegas_mont(n*sizeof(uint64_t));
    CudaBuffer iomegas_mont(n*sizeof(uint64_t));
    // 赋值

    std::vector<uint64_t> tmp(n*batch_size);
    for(size_t i=0; i<tmp.size(); i++) tmp[i] = i;
    buf1.copy_from_host(tmp.data());

    std::vector<uint64_t> omegas_cpu(n);
    omegas_cpu[0] = 1;
    for(size_t i=1; i<n; i++)omegas_cpu[i] = mod_mul(omegas_cpu[i-1], omega, q);
    for(auto& i: omegas_cpu) i = mm.encode(i);  // 转为蒙哥马利
    omegas_mont.copy_from_host(omegas_cpu.data());

    omegas_cpu[0] = 1;
    for(size_t i=1; i<n; i++)omegas_cpu[i] = mod_mul(omegas_cpu[i-1], iomega, q);
    for(auto& i: omegas_cpu) i = mm.encode(i);  // 转为蒙哥马利
    iomegas_mont.copy_from_host(omegas_cpu.data());
    std::cout << "预热与正确性检查" << std::endl;
    float ms1 = cuda_ntt(buf1, omegas_mont, logn, mm, batch_size, true);    // 正向
    float ms2 = cuda_ntt(buf1, iomegas_mont, logn, mm, batch_size, false);  // 逆向

    // 检查结果正确性
    buf1.copy_to_host(tmp.data());
    for(size_t i=0; i<tmp.size(); i++){
        assert(tmp[i] == mod_mul(i, n, q)); // 乘n是因为ntt有乘n的副作用
    }
    std::cout << "正确性检查通过" << std::endl;

    for(int i=1; i<=10; i++)
    {
        std::cout << "第" << i << "次（性能）" << std::endl;
        ms1 = cuda_ntt(buf1, omegas_mont, logn, mm, batch_size, true);
        ms2 = cuda_ntt(buf1, iomegas_mont, logn, mm, batch_size, false);
        std::cout 
            << "(ms1=" << ms1
            << "\tms2=" << ms2
            << "\tsum=" << ms1 + ms2 << ")" << std::endl;
    }
    
    

}
