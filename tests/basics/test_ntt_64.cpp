#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"
#include "ziq_array.hpp"
#include <iostream>


int test_ntt_64()
{
    printf("test_ntt_64\n");
    // n=8 q=1601
    // root-n: 408
    u_int64_t a[8] = {1,2,3,4,5,6,7,8}, b[8], roots[8];
    u_int64_t q=1601, root=408;
    // 需要提供root的各个次幂
    roots[0] = 1;
    for(int i=1;i<8;i++)roots[i] = mod_mul(roots[i-1], root, q);

    ntt_standard_64_with_roots(a, b, roots, 8, q);
    // 预期输出：[36, 1365, 156, 1045, 1597, 548, 1437, 228]
    u_int64_t expected[8] = {36, 1365, 156, 1045, 1597, 548, 1437, 228};

    int first_error = -1;
    for(int i=0;i<8;i++)
    {
        if(b[i] != expected[i])
        {
            std::cout << "error at idx="<< i << ", expect "<<expected[i]<<" but got "<<b[i]<<std::endl;
            first_error = i;
            break;
        }
    }

    if (first_error == -1)
    {
        printf("test_ntt_64 pass!\n");
        return 1;
    }
    else
    {
        printf("test_ntt_64 error...\n");
        return 0;
    }
}

