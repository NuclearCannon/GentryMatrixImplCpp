#include "GentryPoly.hpp"

void GentryPoly::drop_modulus(uint64_t modulus)
{
    assert(is_cpu());
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
    cpu_components().erase(cpu_components().begin() + mod_idx);
    moduli_.erase(moduli_.begin() + mod_idx);
}