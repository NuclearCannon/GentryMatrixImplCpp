#include "u64_array.hpp"

CRTArray CRTArray::mod_reduce(std::shared_ptr<const U64CtxChain> cc2) const
{
    // cc2肯定比cc少一位
    assert(cc2->get_chain_length() == cc_->get_chain_length()-1);
    // cc2和cc的前几位肯定相同
    for(int i=0;i<cc2->get_chain_length();i++)
    {
        assert(cc2->get_mods()[i] == cc_->get_mods()[i]);
    }
    // 可以开始降低模数
    // 首先，减去余数
    vv64 data(data_);
    vec64 remainder = data.back();
    data.pop_back();
    auto& ctxs = cc2->get_ctx();
    auto& mods = cc2->get_mods();

    for(int i=0;i<cc2->get_chain_length(); i++)
    {
        ctxs[i]->sub(data[i], data[i], remainder);
    }
    // 然后，乘以额外模数的乘法逆元
    u64 qo = cc_->get_mods().back();
    for(int i=0;i<cc2->get_chain_length(); i++)
    {
        u64 qo_inv = mod_inv(qo, mods[i]);
        ctxs[i]->mul_scalar(data[i], data[i], qo_inv);
    }
    // 对于余数大于一半的位，+1
    u64 half = qo/2;
    for(int i=0; i<cc2->get_size(); i++)
    {
        if(remainder[i]>half)
        {
            for(int j=0;j<cc2->get_chain_length(); j++)
            {
                data[j][i] += 1;
                data[j][i] %= mods[j];
            }
        }
    }
    return CRTArray(data, cc2);
}