#include "GentryPoly.hpp"

void test_ntt()
{
    size_t n = 256, p = 17;
    uint64_t q = 70368747120641, qr = 6;

    GPComponent x(n, p, q);
    GPCCtx ctx(n, p, q, qr);
    GPComponent y = x;
    y.ntt(ctx);
    y.intt(ctx);
    assert(y.eq(x));
    printf("GPComponent NTT可逆测试成功\n");
}

void test_ntt_cuda()
{
    size_t n = 256, p = 17;
    uint64_t q = 70368747120641, qr = 6;

    GPComponent x(n, p, q), z(n, p, q);
    GPCCtx ctx(n, p, q, qr);
    GPComponentCuda y(n, p, q);
    y.set_from_cpu(x);

    y.ntt(ctx);
    y.intt(ctx);
    y.to_cpu(z);
    assert(z.eq(x));
    printf("GPComponentCuda NTT可逆测试成功\n");
}