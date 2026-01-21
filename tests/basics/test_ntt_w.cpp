#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"
#include "ziq_array.hpp"
#include <iostream>



int test_ntt_w()
{
    printf("test_ntt_w\n");
    // p=5 q=1601 g=2
    // eta: 442
    // [20 741 683 152]
    std::vector<std::string> arr_str = {
        "1", "2", "3", "4"
    };
    fmpz_vector a(arr_str);
    fmpz_t q, eta;
    fmpz_init(q);
    fmpz_init(eta);

    string_to_fmpz("1601", q);
    string_to_fmpz("442", eta);

    TwistedNtterW ntter(
        5, 2, q, eta
    );

    fmpz_vector b(4), c(4);
    ntter.ntt(a, b);    // b=NTT(a)
    ntter.intt(b, c);   // c=iNTT(b)
    std::vector<std::string> result_b = b.export_to_vec_str();
    std::vector<std::string> result_c = c.export_to_vec_str();

    // 预期输出：[20 741 683 152]
    std::vector<std::string> expected = {"20", "741", "683", "152"};

    bool error = 0;
    for(int i=0;i<4;i++)
    {
        if(result_b[i] != expected[i])
        {
            std::cout << "1.error at idx="<< i << ", expect "<<expected[i]<<" but got "<<result_b[i]<<std::endl;
            error = 1;
            break;
        }
    }
    // 预期c==a
    for(int i=0;i<4;i++)
    {
        if(result_c[i] != arr_str[i])
        {
            std::cout << "2.error at idx="<< i << ", expect "<<arr_str[i]<<" but got "<<result_c[i]<<std::endl;
            error = 1;
            break;
        }
    }

    fmpz_clear(q);
    fmpz_clear(eta);

    if (!error)
    {
        printf("test_ntt_w pass!\n");
        return 1;
    }
    else
    {
        printf("test_ntt_w error...\n");
        return 0;
    }
}

