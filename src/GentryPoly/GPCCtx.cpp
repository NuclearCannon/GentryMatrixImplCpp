#include "GentryPoly.hpp"
#include "modops.hpp"


GPCCtx::GPCCtx(size_t n, size_t p, uint64_t q, uint64_t qroot):
    n_(n),
    p_(p),
    q_(q),
    ntter_p_(n, q, qroot),
    ntter_n_(n, q, mod_inv(qroot, q)),
    ntter_w_(p, q, qroot),
    ntter_i_(q, qroot)
{

}