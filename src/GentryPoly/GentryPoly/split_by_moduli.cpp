#include "GentryPoly.hpp"
#include "modops.hpp"
/*
    auto& self_comp = cpu_components();
    const auto& r = self_comp[mod_idx].get_data();
    // 设置remainder
    for(auto& comp : remainder.cpu_components())
    {
        comp.set_from_data(r);
    }
    // 令self-=remainder
    GentryPoly::sub(*this, *this, remainder);
    // 乘以modulus的乘法逆元。除了mod_idx（它应该已经被清零了）
    for(int i=0; i<moduli_.size(); i++)
    {
        if (i == mod_idx)
        {
            // 检查一下是不是清零了
            for(uint64_t j:self_comp[i].get_data())assert(j == 0);
        }
        else
        {
            // 乘以modulus的乘法逆元
            GPComponent::mul_scalar(self_comp[i], self_comp[i], mod_inv(modulus % moduli_[i], moduli_[i]));
        }
    }
*/

std::vector<std::vector<uint64_t>> GentryPoly::split_by_moduli()
{
    assert(is_cpu());
    std::vector<std::vector<uint64_t>> result;
    auto& comp = cpu_components();

    int crt_len = comp.size();


    for(int i=0; i<crt_len; i++)
    {
        GPComponent r = std::move(comp[i]);
        // 给其余位减去r
        for(int j=i+1; j<crt_len; j++)GPComponent::sub(comp[j], comp[j], r);
        // 给其他位乘以乘法逆元
        for(int j=i+1; j<crt_len; j++)GPComponent::mul_scalar(comp[j], comp[j], mod_inv(r.get_q(), comp[j].get_q()));
        // push r
        result.push_back(std::move(r.data_));
    }
    return result;
}