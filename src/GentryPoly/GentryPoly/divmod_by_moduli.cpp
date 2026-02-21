#include "GentryPoly.hpp"
#include "modops.hpp"

void GentryPoly::divmod_by_modulus(
    GentryPoly& remainder,
    uint64_t modulus
)
{
    assert(is_cpu());
    assert(remainder.is_cpu());
    assert(like(remainder));

    // 首先，检查modulus是否确实是自己的一个模数
    int mod_idx = -1;
    for(int i=0; i<moduli_.size(); i++)
    {
        if(moduli_[i] == modulus)
        {
            mod_idx = i;
            break;
        }
    }
    assert(mod_idx != -1);
    // 已经确认modulus确实是自己的一个模数
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
}


void GentryPoly::moduli_extend_mult(uint64_t mod)
{
    assert(is_cpu());
    GentryPoly::mul_scalar(*this, *this, mod);
    this->moduli_.push_back(mod);
    this->cpu_components().push_back(
        GPComponent(n(), p(), mod)
    );
}
void GentryPoly::moduli_extend_unsafe(uint64_t mod)
{
    assert(is_cpu());
    assert(abs() != -1);
    this->moduli_.push_back(mod);
    GPComponent new_comp = GPComponent(n(), p(), mod);
    auto data = cpu_components().front().to_signed();
    this->cpu_components().push_back(
        GPComponent::from_signed_data(n(), p(), mod, data)
    );
}