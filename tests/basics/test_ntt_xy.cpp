#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"
#include "ziq_array.hpp"
#include <iostream>


int test_ntt_xy()
{
    printf("test_ntt_xy\n");
    // n=8 q=1601
    // root-4n: 356
    std::vector<std::string> arr_str = {
        "1", "2", "3", "4", "5", "6", "7", "8"
    };
    fmpz_vector a(arr_str);
    fmpz_t q, zeta;
    fmpz_init(q);
    fmpz_init(zeta);

    string_to_fmpz("1601", q);
    string_to_fmpz("356", zeta);

    TwistedNtterXY ntter(
        8, q, zeta
    );

    fmpz_vector b(8), c(8);

    ntter.ntt(a, b);    // b=NTT(a)
    ntter.intt(b, c);   // c=iNTT(b)


    std::vector<std::string> result_b = b.export_to_vec_str();
    std::vector<std::string> result_c = c.export_to_vec_str();

    // 预期输出：[1100, 1246, 740, 1028, 1189, 1322, 332, 1056]
    std::vector<std::string> expected = {
        "1100", "1246", "740", "1028", "1189", "1322", "332", "1056"
    };

    bool error = 0;
    for(int i=0;i<8;i++)
    {
        if(result_b[i] != expected[i])
        {
            std::cout << "error at idx="<< i << ", expect "<<expected[i]<<" but got "<<result_b[i]<<std::endl;
            error = 1;
            break;
        }
    }
    // 预期c==a
    for(int i=0;i<8;i++)
    {
        if(result_c[i] != arr_str[i])
        {
            std::cout << "error at idx="<< i << ", expect "<<arr_str[i]<<" but got "<<result_c[i]<<std::endl;
            error = 1;
            break;
        }
    }

    fmpz_clear(q);
    fmpz_clear(zeta);

    if (!error)
    {
        printf("test_ntt_xy pass!\n");
        return 1;
    }
    else
    {
        printf("test_ntt_xy error...\n");
        return 0;
    }
}

