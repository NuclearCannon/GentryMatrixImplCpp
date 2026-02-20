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

// ===== Static Operations =====

void GentryPoly::add(GentryPoly& dst, const GentryPoly& a, const GentryPoly& b) {
    if (!a.like(b) || !dst.like(a)) {
        throw std::invalid_argument("GentryPoly::add: operands not like");
    }
    if (a.is_cuda()) {
        auto& d = dst.cuda_components();
        const auto& x = a.cuda_components();
        const auto& y = b.cuda_components();
        for (size_t i = 0; i < d.size(); ++i) {
            GPComponentCuda::add(d[i], x[i], y[i]);
        }
    } else {
        auto& d = dst.cpu_components();
        const auto& x = a.cpu_components();
        const auto& y = b.cpu_components();
        for (size_t i = 0; i < d.size(); ++i) {
            GPComponent::add(d[i], x[i], y[i]);
        }
    }
}

void GentryPoly::neg(GentryPoly& dst, const GentryPoly& a) {
    if (!dst.like(a)) {
        throw std::invalid_argument("GentryPoly::neg: dst not like src");
    }
    if (a.is_cuda()) {
        auto& d = dst.cuda_components();
        const auto& x = a.cuda_components();
        for (size_t i = 0; i < d.size(); ++i) {
            GPComponentCuda::neg(d[i], x[i]);
        }
    } else {
        auto& d = dst.cpu_components();
        const auto& x = a.cpu_components();
        for (size_t i = 0; i < d.size(); ++i) {
            GPComponent::neg(d[i], x[i]);
        }
    }
}

void GentryPoly::mul(GentryPoly& dst, const GentryPoly& a, const GentryPoly& b) {
    if (!a.like(b) || !dst.like(a)) {
        throw std::invalid_argument("GentryPoly::mul: operands not like");
    }
    if (a.is_cuda()) {
        auto& d = dst.cuda_components();
        const auto& x = a.cuda_components();
        const auto& y = b.cuda_components();
        for (size_t i = 0; i < d.size(); ++i) {
            GPComponentCuda::mul(d[i], x[i], y[i]);
        }
    } else {
        auto& d = dst.cpu_components();
        const auto& x = a.cpu_components();
        const auto& y = b.cpu_components();
        for (size_t i = 0; i < d.size(); ++i) {
            GPComponent::mul(d[i], x[i], y[i]);
        }
    }
}

void GentryPoly::ntt(GentryPoly& poly, const GentryPolyCtx& ctx_set) {
    if (poly.is_cuda()) {
        auto& comps = poly.cuda_components();
        const auto& moduli = poly.moduli_;
        for (size_t i = 0; i < comps.size(); ++i) {
            uint64_t q = moduli[i];
            const GPCCtx& ctx = ctx_set.get_ctx(q);
            comps[i].ntt(ctx);
        }
    } else {
        auto& comps = poly.cpu_components();
        const auto& moduli = poly.moduli_;
        for (size_t i = 0; i < comps.size(); ++i) {
            uint64_t q = moduli[i];
            const GPCCtx& ctx = ctx_set.get_ctx(q);
            comps[i].ntt(ctx);
        }
    }
}

void GentryPoly::intt(GentryPoly& poly, const GentryPolyCtx& ctx_set) {
    if (poly.is_cuda()) {
        auto& comps = poly.cuda_components();
        const auto& moduli = poly.moduli_;
        for (size_t i = 0; i < comps.size(); ++i) {
            uint64_t q = moduli[i];
            const GPCCtx& ctx = ctx_set.get_ctx(q);
            comps[i].intt(ctx);
        }
    } else {
        auto& comps = poly.cpu_components();
        const auto& moduli = poly.moduli_;
        for (size_t i = 0; i < comps.size(); ++i) {
            uint64_t q = moduli[i];
            const GPCCtx& ctx = ctx_set.get_ctx(q);
            comps[i].intt(ctx);
        }
    }
}