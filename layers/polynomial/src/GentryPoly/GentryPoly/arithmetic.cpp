#include "GentryPoly.hpp"



template<GentryPoly::CpuOp2 cpu_op, GentryPoly::CudaOp2 cuda_op>
void GentryPoly::_op2_tmpl(GentryPoly& dst, const GentryPoly& a) {
    if (!dst.like(a)) {
        throw std::invalid_argument("GentryPoly: dst not like src");
    }
    if (a.is_cuda()) {
        auto& d = dst.cuda_components();
        const auto& x = a.cuda_components();
        for (size_t i = 0; i < d.size(); ++i) {
            cuda_op(d[i], x[i]);
        }
    } else {
        auto& d = dst.cpu_components();
        const auto& x = a.cpu_components();
        for (size_t i = 0; i < d.size(); ++i) {
            cpu_op(d[i], x[i]);
        }
    }
}


template<GentryPoly::CpuOp3 cpu_op, GentryPoly::CudaOp3 cuda_op>
void GentryPoly::_op3_tmpl(GentryPoly& dst, const GentryPoly& a, const GentryPoly& b) {
    if (!a.like(b) || !dst.like(a)) {
        throw std::invalid_argument("GentryPoly::_op3_tmpl: operands not like");
    }
    if (a.is_cuda()) {
        auto& d = dst.cuda_components();
        const auto& x = a.cuda_components();
        const auto& y = b.cuda_components();
        for (size_t i = 0; i < d.size(); ++i) {
            cuda_op(d[i], x[i], y[i]);
        }
    } else {
        auto& d = dst.cpu_components();
        const auto& x = a.cpu_components();
        const auto& y = b.cpu_components();
        for (size_t i = 0; i < d.size(); ++i) {
            cpu_op(d[i], x[i], y[i]);
        }
    }
}

void GentryPoly::neg(GentryPoly& dst, const GentryPoly& src) {
    _op2_tmpl<GPComponent::neg, GPComponentCuda::neg>(dst, src);
}
void GentryPoly::mont_encode(GentryPoly& dst, const GentryPoly& src) {
    _op2_tmpl<GPComponent::mont_encode, GPComponentCuda::mont_encode>(dst, src);
}
void GentryPoly::mont_decode(GentryPoly& dst, const GentryPoly& src) {
    _op2_tmpl<GPComponent::mont_decode, GPComponentCuda::mont_decode>(dst, src);
}
void GentryPoly::add(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2) {
    _op3_tmpl<GPComponent::add, GPComponentCuda::add>(dst, src1, src2);
}
void GentryPoly::sub(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2) {
    _op3_tmpl<GPComponent::sub, GPComponentCuda::sub>(dst, src1, src2);
}
void GentryPoly::mul(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2) {
    _op3_tmpl<GPComponent::mul, GPComponentCuda::mul>(dst, src1, src2);
}
void GentryPoly::mul_mont(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2) {
    _op3_tmpl<GPComponent::mul_mont, GPComponentCuda::mul_mont>(dst, src1, src2);
}
void GentryPoly::mul_scalar(GentryPoly& dst, const GentryPoly& src1, uint64_t src_scalar) {
    if (!dst.like(src1)) {
        throw std::invalid_argument("GentryPoly: dst not like src");
    }
    if (src1.is_cuda()) {
        auto& d = dst.cuda_components();
        const auto& x = src1.cuda_components();
        for (size_t i = 0; i < d.size(); ++i) {
            GPComponentCuda::mul_scalar(d[i], x[i], src_scalar);
        }
    } else {
        auto& d = dst.cpu_components();
        const auto& x = src1.cpu_components();
        for (size_t i = 0; i < d.size(); ++i) {
            GPComponent::mul_scalar(d[i], x[i], src_scalar);
        }
    }
}

bool GentryPoly::eq(const GentryPoly& other) const
{
    if(other.is_cuda())return eq(other.to_cpu());
    if(is_cuda())return to_cpu().eq(other);

    assert(this->is_cpu());
    assert(other.is_cpu());
    assert(like(other));
    auto& comp1 = cpu_components();
    auto& comp2 = other.cpu_components();
    int len = comp1.size();
    for(int i=0; i<len; i++)
    {
        if(!(comp1[i].eq(comp2[i])))
        {
            return false;
        }
    }
    return true;
}

uint64_t GentryPoly::hash() const
{
    if(is_cuda())return to_cpu().hash();
    uint64_t x = 0;
    for(auto& i: cpu_components())
    {
        for (auto j: i.data_)
        {
            x += j;
        }
    }
    return x;
}