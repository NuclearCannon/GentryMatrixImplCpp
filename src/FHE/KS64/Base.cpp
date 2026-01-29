#include "FHE/key_switch_64.hpp"

KeySwitchKey64Base::KeySwitchKey64Base(const CRTArray& sk_from, const CRTArray& sk_to, u64 B, u64 L, u64 qo, u64 qor):
    B_(B),
    L_(L)
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
    CRTArray sk_to_hig_ntt = sk_to_hig.all_ntt();
    // 计算sk_from * qo
    fmpz_vector sk_from_qo(size);
    _fmpz_vec_scalar_mul_ui(sk_from_qo.raw(), sk_from_fmpz.raw(), size, qo);
    
    for(int l=0; l<L; l++)
    {
        // 加密sk_from_qo
        CRTArray sk_from_qo_Bi = CRTArray::from_fmpz_vector(sk_from_qo, cc_hig_);
        cts_.push_back(encrypt64_CNNN(sk_from_qo_Bi, sk_to_hig_ntt));
        // sk_from_qo *= B
        _fmpz_vec_scalar_mul_ui(sk_from_qo.raw(), sk_from_qo.raw(), size, B);
    }




}
std::pair<CRTArray, CRTArray> KeySwitchKey64Base::key_switch_big_1(const CRTArray& a) const
{
    // 将a转化为fmpz形式
    fmpz_vector a2 = a.to_fmpz_vector();    // 不必center
    
    // 对ct_a_from沿B分解
    int len = a2.len();

    std::vector<CRTArray> a_split;
    fmpz_scalar B(B_);
    while(1)
    {
        // 检查a是不是还有0？
        if (_fmpz_vec_is_zero(a2.raw(), len)) break;
        // 进行一次取模
        fmpz_vector t_mod(len);
        // let t_mod = a % B
        _fmpz_vec_scalar_mod_fmpz(t_mod.raw(), a2.raw(), len, B.raw());
        // let a = a // B
        _fmpz_vec_scalar_tdiv_q_fmpz(a2.raw(), a2.raw(), len, B.raw());
        a_split.push_back(CRTArray::from_fmpz_vector(t_mod, cc_hig_).all_ntt());
    }
    // 检查长度
    if (a_split.size() > L_)
    {
        throw std::invalid_argument("ct_a_from不能被完整分解！\n");
    }
    // 进行一个分别的乘
    CRTArray a_sum = CRTArray::zeros(cc_hig_), b_sum = CRTArray::zeros(cc_hig_);
    for(int i=0;i<a_split.size();i++)
    {
        auto& [cta, ctb] = cts_[i];
        a_sum.adde(a_split[i].mul(cta));
        b_sum.adde(a_split[i].mul(ctb));
    }
    CRTArray a_res = a_sum.all_intt();
    CRTArray b_res = b_sum.all_intt();
    // 降低模数
    return std::make_pair(
        a_res.mod_reduce(cc_low_),
        b_res.mod_reduce(cc_low_)
    );
}


std::pair<CRTArray, CRTArray> KeySwitchKey64Base::key_switch_big_2(const CRTArray& a, const CRTArray& b) const
{
    auto [c, d] = key_switch_big_1(a);
    return std::make_pair(std::move(c), d.add(b));
}