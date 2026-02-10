#include <iostream>
#include <gperftools/profiler.h>

int test_ks64(bool test_base, bool test_crt);
int test_matmul();


int main(int argc, char* argv[])
{
    int sum = 0;
    bool test_base = false;
    bool test_crt = false;
    bool test_mat = false;

    // 解析命令行参数
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-b") {
            test_base = true;
        }
        if (std::string(argv[i]) == "-c") {
            test_crt = true;
        }
        if (std::string(argv[i]) == "-m") {
            test_mat = true;
        }
    }
    if (test_base || test_crt)
    {
        test_ks64(test_base, test_crt);
    }

    if (test_mat)
    {
        test_matmul();
    }
}