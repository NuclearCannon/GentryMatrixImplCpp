#include "GentryPoly.hpp"
#include "matmul.hpp"


void GentryPoly::circledast(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2)
{
    if(dst.is_cpu())
    {
        assert(src1.is_cpu());
        assert(src2.is_cpu());

        auto& comp0 = dst.cpu_components();
        auto& comp1 = src1.cpu_components();
        auto& comp2 = src2.cpu_components();

        for(int i=0; i<comp0.size(); i++)
        {
            circledast_u64_cpu(
                comp0[i].data_.data(), comp1[i].data_.data(), comp2[i].data_.data(), 
                comp0[i].get_n(), comp0[i].get_p(), comp0[i].get_q()
            );
        }
    }
    else
    {
        
        assert(src1.is_cuda());
        assert(src2.is_cuda());
        auto& comp0 = dst.cuda_components();
        auto& comp1 = src1.cuda_components();
        auto& comp2 = src2.cuda_components();

        for(int i=0; i<comp0.size(); i++)
        {
            circledast_u64_gpu(
                comp0[i].data_, comp1[i].data_, comp2[i].data_, 
                comp0[i].get_n(), comp0[i].get_p(), comp0[i].get_q()
            );
        }
        
    }
    
}