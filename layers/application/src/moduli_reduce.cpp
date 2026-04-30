#include "application.hpp"

Ciphertext Ciphertext::moduli_reduce(const std::vector<uint64_t>& moduli) const
{
    GentryPoly a = *a_, b=*b_;
    for(auto i: moduli)
    {
        a.moduli_reduce(i);
        b.moduli_reduce(i);
    }
    
    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a)),
        std::make_unique<GentryPoly>(std::move(b))
    );
}