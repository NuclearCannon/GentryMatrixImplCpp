#include "GentryPoly.hpp"

// GPComponentCuda::

template<GentryPoly::CpuNttOp cpu_op, GentryPoly::CudaNttOp cuda_op>
void GentryPoly::_ntt_tmpl(const GentryPolyCtx&)
{
    if (is_cuda()) {
        auto& comps = cuda_components();
        const auto& moduli = moduli_;
        for (size_t i = 0; i < comps.size(); ++i) {
            uint64_t q = moduli[i];
            const GPCCtx& ctx = ctx_set.get_ctx(q);
            comps[i].cuda_op(ctx);
        }
    } else {
        auto& comps = cpu_components();
        const auto& moduli = moduli_;
        for (size_t i = 0; i < comps.size(); ++i) {
            uint64_t q = moduli[i];
            const GPCCtx& ctx = ctx_set.get_ctx(q);
            comps[i].cpu_op(ctx);
        }
    }
}


void GentryPoly::i_ntt(const GentryPolyCtx& ctx_set) {
    _ntt_tmpl<GPComponent::i_ntt, GPComponentCuda::i_ntt>(ctx_set);
}
void GentryPoly::i_intt(const GentryPolyCtx& ctx_set) {
    _ntt_tmpl<GPComponent::i_intt, GPComponentCuda::i_intt>(ctx_set);
}
void GentryPoly::w_ntt(const GentryPolyCtx& ctx_set) {
    _ntt_tmpl<GPComponent::w_ntt, GPComponentCuda::w_ntt>(ctx_set);
}
void GentryPoly::w_intt(const GentryPolyCtx& ctx_set) {
    _ntt_tmpl<GPComponent::w_intt, GPComponentCuda::w_intt>(ctx_set);
}
void GentryPoly::x_ntt(const GentryPolyCtx& ctx_set) {
    _ntt_tmpl<GPComponent::x_ntt, GPComponentCuda::x_ntt>(ctx_set);
}
void GentryPoly::x_intt(const GentryPolyCtx& ctx_set) {
    _ntt_tmpl<GPComponent::x_intt, GPComponentCuda::x_intt>(ctx_set);
}
void GentryPoly::y_ntt(const GentryPolyCtx& ctx_set) {
    _ntt_tmpl<GPComponent::y_ntt, GPComponentCuda::y_ntt>(ctx_set);
}
void GentryPoly::y_intt(const GentryPolyCtx& ctx_set) {
    _ntt_tmpl<GPComponent::y_intt, GPComponentCuda::y_intt>(ctx_set);
}
