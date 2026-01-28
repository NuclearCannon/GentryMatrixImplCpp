#include <iostream>


int test_enc64();

int main()
{
    int sum = 0;
    sum += test_enc64();

    std::cout << "Total Pass: " << sum << std::endl; 
}