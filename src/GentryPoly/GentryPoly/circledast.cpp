#include "GentryPoly.hpp"
#include "matmul.hpp"


void GentryPoly::circledast(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2)
{
    assert(dst.is_cpu());
    assert(src1.is_cpu());
    assert(src2.is_cpu());

    auto& comp0 = dst.cpu_components();
    auto& comp1 = src1.cpu_components();
    auto& comp2 = src2.cpu_components();

    for(int i=0; i<comp0.size(); i++)
    {
        uint64_t q = comp0[i].get_q();
        fmpz_scalar q_fmpz = fmpz_scalar::from_ui(q);

        MatmulContext mc(comp0[i].get_n(), q_fmpz.raw());
        mc.circledast_u64(comp0[i].data_.data(), comp1[i].data_.data(), comp2[i].data_.data(), comp0[i].get_n(), comp0[i].get_p());
    }
}