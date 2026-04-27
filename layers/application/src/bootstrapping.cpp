#include "application.hpp"

Ciphertext Ciphertext::naive_moduli_extend(uint64_t new_mod) const
{
    GentryPoly a = a_->moduli_extend_fmpz(new_mod);
    GentryPoly b = b_->moduli_extend_fmpz(new_mod);
    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a)),
        std::make_unique<GentryPoly>(std::move(b))
    );
}