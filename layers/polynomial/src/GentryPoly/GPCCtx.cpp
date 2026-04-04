#include "GentryPoly.hpp"
#include "modops.hpp"


GPCCtx::GPCCtx(size_t n, size_t p, uint64_t q, uint64_t qroot):
    n_(n),
    p_(p),
    q_(q),
    ntter_p_(std::make_unique<TwistedNtterXY>(n, q, qroot)),
    ntter_n_(std::make_unique<TwistedNtterXY>(n, q, mod_inv(qroot, q))),
    ntter_w_(std::make_unique<TwistedNtterW>(p, q, qroot)),
    ntter_i_(std::make_unique<NTTerI>(q, qroot))
{

}

GPCCtx::GPCCtx(GPCCtx&&) = default;