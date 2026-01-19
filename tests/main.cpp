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


int test_ntt_standard_with_ctx()
{
    printf("test_ntt_standard_with_ctx\n");
    // n=8 q=1601
    // root-n: 408
    std::vector<std::string> arr_str = {
        "1", "2", "3", "4", "5", "6", "7", "8"
    };
    fmpz_vector a(arr_str);
    fmpz_t q, root;
    fmpz_mod_ctx_t ctx;

    fmpz_init(q);
    fmpz_init(root);

    string_to_fmpz("1601", q);
    fmpz_mod_ctx_init(ctx, q);
    string_to_fmpz("408", root);

    fmpz_vector b(8);

    ntt_standard_flint(a, b, root, 8, ctx);

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
        printf("test_ntt_standard_with_ctx pass!\n");
        return 1;
    }
    else
    {
        printf("test_ntt_standard_with_ctx error...\n");
        return 0;
    }
}



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
int main()
{
    int sum = 0;
    sum += test_ntt_standard();
    sum += test_ntt_standard_with_ctx();
    sum += test_ntt_xy();
    sum += test_ntt_w();
    sum += test_ziq_ctx();

    std::cout << "Total Pass: " << sum << std::endl; 
}