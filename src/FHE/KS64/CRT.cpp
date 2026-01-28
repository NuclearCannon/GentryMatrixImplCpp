#include "FHE/key_switch_64.hpp"

/*

在Base方法中，我们把a分解在基(1, B, B^2, ..., B^{L-1})上
在这里，我们把a分解在基(1, q0, q0*q1, q0*q1*q2, ...)上
我们会先对q0取模得到关于1的分量
然后再对q1取模得到关于q0*q1的分量……

注意，当我们对一个关于模数链(q[0], q[1], ..., q[c-1])的数取模q之后
得到的余数是关于模数q[0]的，没有CRT，甚至可以认为没有模数的

呃，总之，你得到的余数都是没有模数的
你可以生成CRT的cts，然后乘以没模数的余数，对吧？

*/

constexpr bool _PRE_NTT_ = false;


KeySwitchKey64CRT::KeySwitchKey64CRT(const CRTArray& sk_from, const CRTArray& sk_to, u64 qo, u64 qor)
{
    cc_low_ = sk_from.get_cc();
    assert(sk_to.get_cc().get() == cc_low_.get());
    // 构造cc high
    vec64 mods = cc_low_->get_mods();
    vec64 roots = cc_low_->get_roots();
    mods.push_back(qo);
    roots.push_back(qor);
    int n = cc_low_->get_n();
    int p = cc_low_->get_p();
    int size = cc_low_->get_size();
    cc_hig_ = std::make_shared<const U64CtxChain>(n, p, mods, roots);
    // 使用sk_to加密sk_from的各个次幂
    // 先把他们转换为fmpz_vector形式吧
    fmpz_vector sk_from_fmpz = sk_from.to_fmpz_vector_centered();
    fmpz_vector sk_to_fmpz   = sk_to.to_fmpz_vector_centered();
    // 计算sk_to的cc_hig_形式
    CRTArray sk_to_hig =CRTArray::from_fmpz_vector(sk_to_fmpz, cc_hig_);
    // 计算sk_from * qo
    fmpz_vector sk_from_qo(size);
    _fmpz_vec_scalar_mul_ui(sk_from_qo.raw(), sk_from_fmpz.raw(), size, qo);
    // 生成cts
    for(int l=0; l<cc_low_->get_chain_length(); l++)
    {
        // 加密sk_from_qo
        CRTArray sk_from_qo_Bi = CRTArray::from_fmpz_vector(sk_from_qo, cc_hig_);
        auto [a, b] = encrypt64(sk_from_qo_Bi, sk_to_hig);
        if constexpr (_PRE_NTT_)
        {
            cts_.push_back({a.all_ntt(), b.all_ntt()});
        }
        else
        {
            cts_.push_back({a, b});
        }
        
        // sk_from_qo *= mods[l]
        _fmpz_vec_scalar_mul_ui(sk_from_qo.raw(), sk_from_qo.raw(), size, mods[l]);
    }



}
std::pair<CRTArray, CRTArray> KeySwitchKey64CRT::key_switch_big_1(const CRTArray& a) const
{
    auto a_split_raw = a.mod_by_modulo();

    // printf("debug 4\n");
    // 进行一个分别的乘
    std::vector<CRTArray> a_sum, b_sum;
    // printf("debug 5\n");
    for(int i=0;i<a_split_raw.size();i++)
    {
        // printf("debug 6\n");
        auto& [cta, ctb] = cts_[i];
        // printf("debug 7\n");
        auto a_split_i = CRTArray(a_split_raw[i], cc_hig_);
        if constexpr (_PRE_NTT_)
        {
            auto ntted = a_split_i.all_ntt();
            a_sum.push_back(ntted.mul(cta));
            b_sum.push_back(ntted.mul(ctb));
        }
        else
        {
            // 用这个就有问题：
            // a_sum.push_back(a_split_i.all_ntt().mul(cta.all_ntt()).all_intt());
            // b_sum.push_back(a_split_i.all_ntt().mul(ctb.all_ntt()).all_intt());
            // 用这个就没问题：
            a_sum.push_back(a_split_i.mul_poly(cta));
            b_sum.push_back(a_split_i.mul_poly(ctb));
        }
        
    }
    if (_PRE_NTT_)
    {
        CRTArray a_res = CRTArray::sum(a_sum).all_intt();
        CRTArray b_res = CRTArray::sum(b_sum).all_intt();
        // 降低模数
        return std::make_pair(
            a_res.mod_reduce(cc_low_),
            b_res.mod_reduce(cc_low_)
        );
    }
    else
    {
        CRTArray a_res = CRTArray::sum(a_sum);
        CRTArray b_res = CRTArray::sum(b_sum);
        // 降低模数
        return std::make_pair(
            a_res.mod_reduce(cc_low_),
            b_res.mod_reduce(cc_low_)
        );
    }

}


std::pair<CRTArray, CRTArray> KeySwitchKey64CRT::key_switch_big_2(const CRTArray& a, const CRTArray& b) const
{
    auto [c, d] = key_switch_big_1(a);
    return std::make_pair(std::move(c), d.add(b));
}