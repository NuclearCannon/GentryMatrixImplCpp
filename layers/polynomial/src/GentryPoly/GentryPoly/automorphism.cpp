#include "GentryPoly.hpp"
#include "modops.hpp"

GentryPoly GentryPoly::automorphism(int x, int y, int w) const
{
    assert(x!=0);
    assert(is_cpu());
    if(x == 1 && y == 1 && w==0)
    {
        return *this;
    }
    if(x == 1 && y == 1)
    {
        GentryPoly res = zeros_like(*this, GPDevice::CPU);
        auto& src = this->cpu_components();
        auto& dst = res.cpu_components();
        for(int i=0; i<src.size(); i++)
        {
            src[i].automorphism_W(dst[i], w);
        }
        return res;
    }
    else if (w == 0)
    {
        GentryPoly res = zeros_like(*this, GPDevice::CPU);
        auto& src = this->cpu_components();
        auto& dst = res.cpu_components();
        for(int i=0; i<src.size(); i++)
        {
            src[i].automorphism_XY(dst[i], x, y);
        }
        return res;
    }
    else
    {
        return automorphism(x, y, 0).automorphism(1, 1, w);
    }
    
}

GentryPoly GentryPoly::rotate(int x, int y, int w) const
{
    int n = this->n();
    return automorphism(mod_pow(5, x, 4*n), mod_pow(5, y, 4*n), w);
}