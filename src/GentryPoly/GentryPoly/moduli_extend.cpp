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