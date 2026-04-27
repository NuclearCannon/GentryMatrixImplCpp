#include "application.hpp"

Ciphertext Ciphertext::moduli_reduce(uint64_t mod) const
{
    GentryPoly a = *a_, b=*b_;
    a.moduli_reduce(mod);
    b.moduli_reduce(mod);
    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a)),
        std::make_unique<GentryPoly>(std::move(b))
    );
}