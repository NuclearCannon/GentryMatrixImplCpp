#include "flints.hpp"

std::string fmpz_to_string(fmpz_t x)
{
    char* str = fmpz_get_str(nullptr, 10, x); // 获取十进制字符串
    std::string result(str);
    flint_free(str); // 释放内存
    return result;
}


void string_to_fmpz(const std::string s, fmpz_t x)
{
    fmpz_set_str(x, s.c_str(), 10);
}
