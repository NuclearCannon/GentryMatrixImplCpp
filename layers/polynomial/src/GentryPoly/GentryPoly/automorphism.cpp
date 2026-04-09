#include "GentryPoly.hpp"

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