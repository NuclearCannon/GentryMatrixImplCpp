#include "u64_array.hpp"

vv64 CRTArray::mod_by_modulo() const
{
    int size = cc_->get_size();
    int chain_len = cc_->get_chain_length();
    auto& ctxs = cc_->get_ctx();
    auto& mods = cc_->get_mods();

    vv64 result(chain_len);
    vv64 data = data_;

    for(int i=0; i<chain_len; i++)
    {
        // 现在开始分离：原式 % (prod_{j=0}^{i} q[j])
        vec64& remainder = data[i];
        // 令剩余部分减去remainder
        for(int j=i+1; j<chain_len; j++)
        {
            ctxs[j]->sub(data[j], data[j], remainder);
        }
        // 令剩余部分除以mod[i]
        for(int j=i+1; j<chain_len; j++)
        {
            uint64_t inv = mod_inv(mods[i], mods[j]);
            ctxs[j]->mul_scalar(data[j], data[j], inv);
        }
        // push
        result[i].swap(remainder);
    }
    return result;
}

