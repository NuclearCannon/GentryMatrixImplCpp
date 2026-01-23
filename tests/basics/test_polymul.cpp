#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"
#include "ziq_array.hpp"
#include <iostream>


int test_polymul()
{
    printf("test_polymul\n");
    // n=8 p=5 q=1601
    // zeta: 356  eta: 442
    int n=8, p=5, g=2;
    int size_ = 2*(p-1)*n*n;
    fmpz_scalar q("1601");
    fmpz_scalar root_q("3");
    ZiqArrayContext ctx(n, p, g, q.raw(), root_q.raw());
    fmpz_vector data(size_);
    for(int i=0;i<size_;i++)fmpz_set_ui(data[i], i);

    ZiqArray A1 = ctx.uniform();
    ZiqArray A2 = ctx.uniform();

    ZiqArray A1_ntt = A1.iw_ntt().xy_ntt();
    ZiqArray A2_ntt = A2.iw_ntt().xy_ntt();
    ZiqArray B_ntt = A1_ntt.mul(A2_ntt);
    ZiqArray B = B_ntt.xy_intt().iw_intt();

    // 朴素计算C
    fmpz_vector c_data(size_);
    fmpz_scalar coeff1, coeff2, coeff3;

    fmpz_mod_ctx_t qctx;
    fmpz_mod_ctx_init(qctx, q.raw());


    for(int i1=0;i1<2;i1++)
    {
        for(int w1=0;w1<p-1;w1++)
        {
            for(int x1=0;x1<n;x1++)
            {
                for(int y1=0;y1<n;y1++)
                {
                    for(int i2=0;i2<2;i2++)
                    {
                        for(int w2=0;w2<p-1;w2++)
                        {
                            for(int x2=0;x2<n;x2++)
                            {
                                for(int y2=0;y2<n;y2++)
                                {
                                    
                                    int i3, w3, x3, y3, sign = +1;
                                    // 分别表示结果中i, W, X, Y的次数
                                    i3 = i1 + i2;
                                    w3 = w1 + w2;
                                    x3 = x1 + x2;
                                    y3 = y1 + y2;
                                    // X^n = i
                                    i3 += x3 / n;
                                    x3 %= n;
                                    // Y^n = -i
                                    i3 -= y3 / n;
                                    y3 %= n;
                                    // W^p=1
                                    w3 %= p;
                                    // i^4 == 1
                                    i3 = (i3+4)%4;
                                    if (i3 >= 2)
                                    {
                                        i3 -= 2;
                                        sign = -1;
                                    }
                                    // 一个新的项目
                                    // 读取系数
                                    fmpz_set(coeff1.raw(), A1.data()[i1*(p-1)*n*n + w1*n*n + x1*n + y1]);
                                    
                                    fmpz_set(coeff2.raw(), A2.data()[i2*(p-1)*n*n + w2*n*n + x2*n + y2]);
                                    // 计算
                                    fmpz_mod_mul(coeff3.raw(), coeff1.raw(), coeff2.raw(), qctx);
                                    if (sign == -1)
                                    {
                                        fmpz_mod_neg(coeff3.raw(), coeff3.raw(), qctx);
                                    }
                                    if (w3 == p-1)
                                    {
                                        for(int k=0;k<p-1;k++)
                                        {
                                            fmpz* dst = c_data[i3*(p-1)*n*n + k*n*n + x3*n + y3];
                                            fmpz_mod_sub(dst, dst, coeff3.raw(), qctx);
                                        }
                                    }
                                    else
                                    {
                                        fmpz* dst = c_data[i3*(p-1)*n*n + w3*n*n + x3*n + y3];
                                        fmpz_mod_add(dst, dst, coeff3.raw(), qctx);
                                    }

                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fmpz_mod_ctx_clear(qctx);
    ZiqArray C(c_data, &ctx);
    ZiqArray E = B.sub(C);

    // 比较朴素乘法结果和NTT乘法结果
    int e = E.data().mod_centered(q.raw()).max_abs();

    if (e == 0)
    {
        printf("test_polymul pass!\n");
        return 1;
    }
    else
    {
        printf("test_polymul error(%d)\n", e);
        return 0;
    }
}

