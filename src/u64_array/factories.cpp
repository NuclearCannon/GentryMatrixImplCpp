#include "u64_array.hpp"
#include "random.hpp"

CRTArray CRTArray::zeros(std::shared_ptr<const U64CtxChain> cc)
{
    return CRTArray(cc);    // 普通构造函数就是全0的
}
CRTArray CRTArray::uniform(std::shared_ptr<const U64CtxChain> cc)
{
    // 构造合适的data
    // 对于mods[i], data[i]是长度为size的[0, mods[i])均匀分布
    vv64 data;
    const vec64& mods = cc->get_mods();
    int size = cc->get_size();
    for(u64 mod : mods)
    {
        vec64 data_i(size);
        for(int j=0;j<size;j++)data_i[j] = ramdom_generators::randu64(0, mod-1);
        data.push_back(std::move(data_i));
    }
    return CRTArray(data, cc);
}
CRTArray CRTArray::dg(std::shared_ptr<const U64CtxChain> cc)
{
    const vec64& mods = cc->get_mods();
    int size = cc->get_size();
    vv64 data(cc->get_chain_length(), vec64(size, 0));  // 先造个空的出来
    for(int i=0;i<size;i++)
    {
        // 生成第i个随机数
        i64 t = ramdom_generators::dg(5);
        // 取余
        for(int j=0;j<mods.size();j++)
        {
            u64 mod = mods[j];
            i64 t2 = t % (i64)mod;
            while(t2<0)t2+=mod;
            u64 t3 = t2;
            data[j][i] = t3 % mod;
        }
    }
    return CRTArray(data, cc);
}

CRTArray CRTArray::sk(std::shared_ptr<const U64CtxChain> cc)
{
    const vec64& mods = cc->get_mods();
    int size = cc->get_size();
    int n = cc->get_n();
    vv64 data(cc->get_chain_length(), vec64(size, 0));  // 先造个空的出来
    for(int i=0;i<size;i+=n)
    {
        // 生成第i个随机数
        i64 t = ramdom_generators::dg(5);
        // 取余
        for(int j=0;j<mods.size();j++)
        {
            u64 mod = mods[j];
            i64 t2 = t % (i64)mod;
            while(t2<0)t2+=mod;
            u64 t3 = t2;
            data[j][i] = t3 % mod;
        }
    }
    return CRTArray(data, cc);
}


CRTArray CRTArray::randint(std::shared_ptr<const U64CtxChain> cc, i64 start, i64 end)
{
    
    const vec64& mods = cc->get_mods();
    int size = cc->get_size();
    vv64 data(cc->get_chain_length(), vec64(size, 0));  // 先造个空的出来
    for(int i=0;i<size;i++)
    {
        // 生成第i个随机数
        i64 t = ramdom_generators::randi64(start, end);
        // 取余
        for(int j=0;j<mods.size();j++)
        {
            u64 mod = mods[j];
            i64 t2 = t % (i64)mod;
            while(t2<0)t2+=mod;
            u64 t3 = t2;
            data[j][i] = t3 % mod;
        }
    }
    return CRTArray(data, cc);
}

