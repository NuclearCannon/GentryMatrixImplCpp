#include <iostream>
#include <gperftools/profiler.h>

int test_ks64();

int main()
{
    int sum = 0;
    ProfilerStart("bench64.prof");
    sum += test_ks64();
    ProfilerStop();
    std::cout << "Total Pass: " << sum << std::endl; 
}