#include <iostream>
#include <gperftools/profiler.h>

int test_ks64(bool test_base, bool test_crt);

int main(int argc, char* argv[])
{
    int sum = 0;
    bool test_base = false;
    bool test_crt = false;

    // 解析命令行参数
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-b") {
            test_base = true;
        }
        if (std::string(argv[i]) == "-c") {
            test_crt = true;
        }
    }

    sum += test_ks64(test_base, test_crt);
    std::cout << "Total Pass: " << sum << std::endl; 
}