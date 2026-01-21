#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"
#include "ziq_array.hpp"
#include <iostream>


int test_ntt_standard()
{
    printf("test_ntt_standard\n");
    // n=8 q=1601
    // root-n: 408
    std::vector<std::string> arr_str = {
        "1", "2", "3", "4", "5", "6", "7", "8"
    };
    fmpz_vector a(arr_str);
    fmpz_t q, root;
    fmpz_init(q);
    fmpz_init(root);

    string_to_fmpz("1601", q);
    string_to_fmpz("408", root);

    fmpz_vector b(8);

    ntt_standard_flint(a, b, root, 8, q);

    std::vector<std::string> result = b.export_to_vec_str();

    // 预期输出：[36, 1365, 156, 1045, 1597, 548, 1437, 228]
    std::vector<std::string> expected = {
        "36", "1365", "156", "1045", "1597", "548", "1437", "228"
    };

    int first_error = -1;
    for(int i=0;i<8;i++)
    {
        if(result[i] != expected[i])
        {
            std::cout << "error at idx="<< i << ", expect "<<expected[i]<<" but got "<<result[i]<<std::endl;
            first_error = i;
            break;
        }
    }

    fmpz_clear(q);
    fmpz_clear(root);

    if (first_error == -1)
    {
        printf("test_ntt_standard pass!\n");
        return 1;
    }
    else
    {
        printf("test_ntt_standard error...\n");
        return 0;
    }
}

