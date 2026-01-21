#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"
#include "ziq_array.hpp"
#include <iostream>


int test_ntt_standard();
int test_ntt_xy();
int test_ntt_w();
int test_ziq_ctx();
int test_circledast();



int main()
{
    int sum = 0;
    sum += test_ntt_standard();
    sum += test_ntt_xy();
    sum += test_ntt_w();
    sum += test_ziq_ctx();
    sum += test_circledast();

    std::cout << "Total Pass: " << sum << std::endl; 
}