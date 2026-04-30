#include "application.hpp"

Ciphertext Ciphertext::naive_moduli_extend(const std::vector<uint64_t>& extra_moduli) const
{
    GentryPoly a = a_->moduli_extend_fmpz(extra_moduli);
    GentryPoly b = b_->moduli_extend_fmpz(extra_moduli);
    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a)),
        std::make_unique<GentryPoly>(std::move(b))
    );
}