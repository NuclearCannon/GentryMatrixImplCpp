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

KeySwitchKey64CRT::KeySwitchKey64CRT(std::shared_ptr<const U64CtxChain> cc_low,std::shared_ptr<const U64CtxChain> cc_hig, std::vector<std::pair<CRTArray, CRTArray>> cts):
    cc_low_(cc_low),
    cc_hig_(cc_hig),
    cts_(std::move(cts))
{

}

KeySwitchKey64CRT KeySwitchKey64CRT::ksk_gen(const CRTArray& sk_from, const CRTArray& sk_to, u64 qo, u64 qor)
{
    std::shared_ptr<const U64CtxChain> cc_low = sk_from.get_cc();
    assert(sk_to.get_cc().get() == cc_low.get());
    // 构造cc high
    vec64 mods = cc_low->get_mods();
    vec64 roots = cc_low->get_roots();
    mods.push_back(qo);
    roots.push_back(qor);
    int n = cc_low->get_n();
    int p = cc_low->get_p();
    int size = cc_low->get_size();
    std::shared_ptr<const U64CtxChain> cc_hig = std::make_shared<const U64CtxChain>(n, p, mods, roots);
    // 使用sk_to加密sk_from的各个次幂
    // 先把他们转换为fmpz_vector形式吧
    fmpz_vector sk_from_fmpz = sk_from.to_fmpz_vector_centered();
    fmpz_vector sk_to_fmpz   = sk_to.to_fmpz_vector_centered();
    // 计算sk_to的cc_hig_形式
    CRTArray sk_to_hig =CRTArray::from_fmpz_vector(sk_to_fmpz, cc_hig);
    CRTArray sk_to_hig_ntt = sk_to_hig.all_ntt();
    // 计算sk_from * qo
    CRTArray sk_from_qo_Bi = CRTArray::from_fmpz_vector(sk_from_fmpz, cc_hig);
    sk_from_qo_Bi.mul_scalar_e(qo);
    // 生成cts
    std::vector<std::pair<CRTArray, CRTArray>> cts;
    for(int l=0; l<cc_low->get_chain_length(); l++)
    {
        // 加密sk_from_qo
        cts.push_back(encrypt64_CNNN(sk_from_qo_Bi, sk_to_hig_ntt));
        // sk_from_qo *= mods[l]
        sk_from_qo_Bi.mul_scalar_e(mods[l]);
    }
    for(auto &p: cts)
    {
        p.first.mont_encode_inplace();
        p.second.mont_encode_inplace();
    }
    return KeySwitchKey64CRT(
        cc_low, cc_hig, std::move(cts)
    );
}
std::pair<CRTArray, CRTArray> KeySwitchKey64CRT::key_switch_big_1(const CRTArray& a) const
{
    auto a_split_raw = a.mod_by_modulo();

    // printf("debug 4\n");
    // 进行一个分别的乘
    CRTArray a_sum = CRTArray::zeros(cc_hig_), b_sum = CRTArray::zeros(cc_hig_);
    CRTArray mul_result = CRTArray::zeros(cc_hig_);
    CRTArray ntted = CRTArray::zeros(cc_hig_);
    CRTArray a_split_i = CRTArray::zeros(cc_hig_);
    CRTArray a_res = CRTArray::zeros(cc_hig_);
    CRTArray b_res = CRTArray::zeros(cc_hig_);
    // printf("debug 5\n");
    for(int i=0;i<a_split_raw.size();i++)
    {
        // printf("debug 6\n");
        auto& [cta, ctb] = cts_[i];
        // printf("debug 7\n");
        a_split_i.set_from_raw(a_split_raw[i]);
        a_split_i.all_ntt_to(ntted);
        // ntted.mont_encode_inplace(); // 理论上这一行应该存在，但是和下面抵消了
        CRTArray::mul_mont3(mul_result, ntted, cta);
        a_sum.adde(mul_result);
        CRTArray::mul_mont3(mul_result, ntted, ctb);
        b_sum.adde(mul_result);
    }
    a_sum.all_intt_to(a_res);
    b_sum.all_intt_to(b_res);
    // 理论上这两行应该存在，但是和上面抵消了
    // a_res.mont_decode_inplace();
    // b_res.mont_decode_inplace();

    // 降低模数
    return std::make_pair(
        a_res.mod_reduce(cc_low_),
        b_res.mod_reduce(cc_low_)
    );


}


std::pair<CRTArray, CRTArray> KeySwitchKey64CRT::key_switch_big_2(const CRTArray& a, const CRTArray& b) const
{
    auto [c, d] = key_switch_big_1(a);
    return std::make_pair(std::move(c), d.add(b));
}