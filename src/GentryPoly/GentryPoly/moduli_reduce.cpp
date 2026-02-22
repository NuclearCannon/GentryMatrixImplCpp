#include "GentryPoly.hpp"
#include "modops.hpp"

void GentryPoly::moduli_reduce(
    uint64_t modulus
)
{
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
    if(is_cpu())
    {
        auto& self_comp = cpu_components();
        auto& r = self_comp[mod_idx];
        for(int i=0; i<moduli_.size(); i++)
        {
            if (i != mod_idx)
            {
                auto& comp = self_comp[i];
                // 乘以modulus的乘法逆元
                GPComponent::sub(comp, comp, r);
                GPComponent::mul_scalar(comp, comp, mod_inv(modulus % moduli_[i], moduli_[i]));
            }
        }
        self_comp.erase(self_comp.begin() + mod_idx);
        moduli_.erase(moduli_.begin() + mod_idx);
    }
    else if (is_cuda())
    {
        auto& self_comp = cuda_components();
        auto& r = self_comp[mod_idx];
        for(int i=0; i<moduli_.size(); i++)
        {
            if (i != mod_idx)
            {
                auto& comp = self_comp[i];
                // 乘以modulus的乘法逆元
                GPComponentCuda::sub(comp, comp, r);
                GPComponentCuda::mul_scalar(comp, comp, mod_inv(modulus % moduli_[i], moduli_[i]));
            }

        }
        self_comp.erase(self_comp.begin() + mod_idx);
        moduli_.erase(moduli_.begin() + mod_idx);
    }
    else
    {
        throw std::runtime_error("moduli_reduce: You should not arrive here\n");
    }

}


