#include "FHE/encrypt.hpp"
#include "ziq_array.hpp"
#include <iostream>


int test_enc();
int test_ks();
int test_circledast_ct();

int main()
{
    int sum = 0;
    // sum += test_enc();
    // sum += test_ks();
    sum += test_circledast_ct();

    std::cout << "Total Pass: " << sum << std::endl; 
}