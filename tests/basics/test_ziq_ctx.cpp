#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"
#include "ziq_array.hpp"
#include <iostream>


int test_ziq_ctx()
{
    printf("test_ziq_ctx\n");
    // n=8 p=5 q=1601
    // zeta: 356  eta: 442
    int n=8, p=5, g=2;
    int size_ = 2*(p-1)*n*n;
    fmpz_t q, zeta, eta;
    fmpz_init(q);
    fmpz_init(zeta);
    fmpz_init(eta);

    string_to_fmpz("1601", q);
    string_to_fmpz("356", zeta);
    string_to_fmpz("442", eta);
    ZiqArrayContext ctx(n, p, g, q, zeta, eta);
    fmpz_vector data(size_);
    for(int i=0;i<size_;i++)fmpz_set_ui(data[i], i);

    int error0 = -1, error1 = -1;

    auto r1 = ctx.iw_ntt(data);
    auto r2 = ctx.iw_intt(r1);
    for(int i=0;i<size_;i++)
    {
        if(fmpz_cmp(r2[i], data[i]) != 0)
        {
            error0 = i;
            break;
        }
    }

    auto r3 = ctx.xy_ntt(data);
    auto r4 = ctx.xy_intt(r3);
    for(int i=0;i<size_;i++)
    {
        if(fmpz_cmp(r4[i], data[i]) != 0)
        {
            error1 = i;
            break;
        }
    }


    fmpz_clear(q);
    fmpz_clear(zeta);
    fmpz_clear(eta);

    if ((error0 == -1) && (error1 == -1))
    {
        printf("test_ziq_ctx pass!\n");
        return 1;
    }
    else
    {
        printf("test_ziq_ctx error(%d, %d)\n", error0, error1);
        return 0;
    }
}

