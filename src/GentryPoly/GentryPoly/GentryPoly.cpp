#include "GentryPoly.hpp"
#include <algorithm>
#include <cassert>

// ======================
// GentryPolyCtx
// ======================

GentryPolyCtx::GentryPolyCtx(size_t n, size_t p, const std::vector<QRootPair>& q_and_qroots)
    : n_(n), p_(p) {
    for (const auto& [q, qroot] : q_and_qroots) {
        ctx_map_.emplace(q, GPCCtx(n, p, q, qroot));
    }
}

const GPCCtx& GentryPolyCtx::get_ctx(uint64_t q) const {
    auto it = ctx_map_.find(q);
    if (it == ctx_map_.end()) {
        throw std::invalid_argument("Modulus q not found in GentryPolyCtx");
    }
    return it->second;
}

// ======================
// GentryPoly
// ======================

GentryPoly::GentryPoly(bool is_cuda, std::vector<uint64_t> moduli,
                       std::vector<GPComponent> cpu_comps)
    : is_cuda_(false), moduli_(std::move(moduli)) {
    storage_ = std::move(cpu_comps);
}

GentryPoly::GentryPoly(bool is_cuda, std::vector<uint64_t> moduli,
                       std::vector<GPComponentCuda> cuda_comps)
    : is_cuda_(true), moduli_(std::move(moduli)) {
    storage_ = std::move(cuda_comps);
}

GentryPoly GentryPoly::from_coeffs(
    size_t n, size_t p,
    const std::vector<uint64_t>& moduli,
    const std::vector<std::vector<uint64_t>>& coeffs_mod_q) {

    if (moduli.size() != coeffs_mod_q.size()) {
        throw std::invalid_argument("moduli and coeffs_mod_q size mismatch");
    }

    std::vector<GPComponent> comps;
    comps.reserve(moduli.size());
    for (size_t i = 0; i < moduli.size(); ++i) {
        uint64_t q = moduli[i];
        const auto& coeffs = coeffs_mod_q[i];
        comps.emplace_back(GPComponent::from_data(n, p, q, coeffs));
    }
    return GentryPoly(false, moduli, std::move(comps));
}

void GentryPoly::set_from_vec64(const std::vector<uint64_t>& v)
{
    assert(is_cpu());
    for(auto& i: cpu_components())
    {
        i.set_from_data(v);
    }
}

size_t GentryPoly::n() const {
    if (moduli_.empty()) return 0;
    if (is_cuda_) return cuda_components()[0].get_n();
    else          return cpu_components()[0].get_n();
}

size_t GentryPoly::p() const {
    if (moduli_.empty()) return 0;
    if (is_cuda_) return cuda_components()[0].get_p();
    else          return cpu_components()[0].get_p();
}

bool GentryPoly::like_no_device(const GentryPoly& other) const {
    if (moduli_ != other.moduli_) return false;
    if (moduli_.empty()) return true;
    // Check n and p via first component
    return n() == other.n() && p() == other.p();
}

bool GentryPoly::like(const GentryPoly& other) const {
    return is_cuda_ == other.is_cuda_ && like_no_device(other);
}

// Accessors
const GentryPoly::CpuStorage& GentryPoly::cpu_components() const {
    return std::get<CpuStorage>(storage_);
}
GentryPoly::CpuStorage& GentryPoly::cpu_components() {
    return std::get<CpuStorage>(storage_);
}
const GentryPoly::CudaStorage& GentryPoly::cuda_components() const {
    return std::get<CudaStorage>(storage_);
}
GentryPoly::CudaStorage& GentryPoly::cuda_components() {
    return std::get<CudaStorage>(storage_);
}


int64_t GentryPoly::abs() const
{
    if(is_cuda())return this->to_cpu().abs();
    const std::vector<GPComponent>& comp = cpu_components();
    size_t size = 2*(p()-1)*n()*n();
    int64_t result = 0;

    for(size_t i=0; i<size; i++)
    {
        int64_t value = 0;
        bool valid = false;

        for(auto& c: comp)
        {
            int64_t ci = c.get_data()[i];
            uint64_t qi = c.get_q();
            assert(ci >= 0);
            assert(ci < qi);
            if(ci*2>qi)ci -= qi;
            if(!valid)
            {
                valid = true;
                value = ci;
            }
            if(value != ci)
            {
                return -1;
            }
        }
        value = std::abs(value);
        if(value>result)result = value;
    }
    return result;
}