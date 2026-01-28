#include <iostream>


int test_enc64();
int test_ks64B();
int test_ks64C();

int main()
{
    int sum = 0;
    // sum += test_enc64();
    // sum += test_ks64B();
    sum += test_ks64C();
    std::cout << "Total Pass: " << sum << std::endl; 
}