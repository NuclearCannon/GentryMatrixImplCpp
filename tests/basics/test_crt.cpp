#include "CRT.hpp"

// 选取几个约为2^60的质数
// 质数: 1152921504606847009       原根: 13
// 质数: 1152921504606847067       原根: 6
// 质数: 1152921504606847081       原根: 11
// 质数: 1152921504606847123       原根: 2
// 质数: 1152921504606847127       原根: 5

// 大质数: 1962189608304824093585857412065247859367422524039

int test_crt()
{
    printf("test_crt\n");
    fmpz_scalar Q("1962189608304824093585857412065247859367422524039");
    std::vector<u64> mods = {
        1152921504606847009UL, 
        1152921504606847067UL, 
        1152921504606847081UL, 
        1152921504606847123UL, 
        1152921504606847127UL
    };

    // 随机生成一段数据
    int len = 32;
    fmpz_vector data = fmpz_vector::uniform(len, Q.raw());
    auto crt_res = crt(data, mods);
    fmpz_vector icrt_res(len);
    icrt(icrt_res, crt_res, mods);
    // 检查复原是否成功
    bool error = false;
    for(int i=0;i<len;i++)
    {
        if(fmpz_cmp(data[i], icrt_res[i]))
        {
            error = true;
            break;
        }
    }
    if (error)
    {
        printf("icrt失败");
        printf("===============================\n");
        data.print();
        printf("===============================\n");
        icrt_res.print();
        printf("===============================\n");
        return 0;
    }
    else
    {
        printf("test_crt pass!\n");
        return 1;
    }
}