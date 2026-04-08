#include "GentryPoly.hpp"

GentryPoly GentryPoly::transpose() const
{
    assert(is_cpu());
    GentryPoly res = zeros_like(*this, GPDevice::CPU);
    auto& src = this->cpu_components();
    auto& dst = res.cpu_components();
    for(int i=0; i<src.size(); i++)
    {
        src[i].transpose(dst[i]);
    }
    return res;
}
GentryPoly GentryPoly::conj() const
{
    assert(is_cpu());
    GentryPoly res = zeros_like(*this, GPDevice::CPU);
    auto& src = this->cpu_components();
    auto& dst = res.cpu_components();
    for(int i=0; i<src.size(); i++)
    {
        src[i].conj(dst[i]);
    }
    return res;
}
GentryPoly GentryPoly::w_inv() const
{
    assert(is_cpu());
    GentryPoly res = zeros_like(*this, GPDevice::CPU);
    auto& src = this->cpu_components();
    auto& dst = res.cpu_components();
    for(int i=0; i<src.size(); i++)
    {
        src[i].w_inv(dst[i]);
    }
    return res;
}