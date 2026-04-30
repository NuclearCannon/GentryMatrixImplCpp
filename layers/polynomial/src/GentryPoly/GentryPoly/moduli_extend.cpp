#include "GentryPoly.hpp"
#include "modops.hpp"


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

GentryPoly GentryPoly::moduli_extend_fmpz(const std::vector<uint64_t>& extra_moduli) const
{
    assert(is_cpu());
    auto new_moduli = moduli_;
    for(auto mod: extra_moduli)new_moduli.push_back(mod);
    auto vec = to_fmpz_vector();
    return GentryPoly::from_fmpz(n(), p(), new_moduli, vec);
}

GentryPoly GentryPoly::moduli_extend_garner(uint64_t mod) const
{
    throw std::runtime_error("moduli_extend_garner 未被实现，请使用 moduli_extend_fmpz 替代\n");
}