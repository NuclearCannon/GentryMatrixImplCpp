#include <iostream>


int test_enc64();
int test_ks64B();


int main()
{
    int sum = 0;
    // sum += test_enc64();
    sum += test_ks64B();;
    std::cout << "Total Pass: " << sum << std::endl; 
}