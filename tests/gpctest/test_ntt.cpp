#include "GentryPoly.hpp"

void test_ntt()
{
    size_t n = 256, p = 17;
    uint64_t q = 70368747120641, qr = 6;

    GPComponent x(256, 17, 70368747120641);
    GPCCtx ctx(256, 17, 70368747120641, 6);
    GPComponent y = x;
    y.ntt(ctx);
    y.intt(ctx);
    assert(y.eq(x));
    printf("GPComponent NTT可逆测试成功\n");
}