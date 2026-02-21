#include "GentryPoly.hpp"
#include "modops.hpp"

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