#include "GentryPoly.hpp"
#include "modops.hpp"

GentryPoly GentryPoly::automorphism2(int x, int y) const
{
    assert(is_cpu());
    GentryPoly res = zeros_like(*this, GPDevice::CPU);
    auto& src = this->cpu_components();
    auto& dst = res.cpu_components();
    for(int i=0; i<src.size(); i++)
    {
        src[i].automorphism_XY(dst[i], x, y);
    }
    return res;
}

GentryPoly GentryPoly::rotate2(int x, int y) const
{
    int n = this->n();
    return automorphism2(mod_pow(5, x, 4*n), mod_pow(5, y, 4*n));
}