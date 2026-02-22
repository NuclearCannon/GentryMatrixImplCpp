#include "FHE/key_switch_gp.hpp"




KeySwitchKeyGP::KeySwitchKeyGP(
    std::vector<std::pair<GentryPoly, GentryPoly>> cts,
    std::vector<std::pair<GentryPoly, GentryPoly>> cts_cuda,
    uint64_t qo
):
    cts_(std::move(cts)), 
    cts_cuda_(std::move(cts_cuda)),
    qo_(qo)
{

}

KeySwitchKeyGP KeySwitchKeyGP::ksk_gen(const GentryPoly& sk_from, const GentryPoly& sk_to, uint64_t qo, const GentryPolyCtx& ctx)
{
    // 这里假设sk_from, sk_to都是含qo模数的，且sk_from已经乘过qo
    auto mods = sk_to.moduli();
    GentryPoly sk_from_cp = sk_from.to_cpu();
    GentryPoly sk_to_cpu = sk_to.to_cpu();
    std::vector<std::pair<GentryPoly, GentryPoly>> cts;
    for(auto mod : mods)
    {
        if(mod == qo)continue;
        cts.push_back(encrypt_gp(sk_from_cp, sk_to_cpu, ctx));
        GentryPoly::mul_scalar(sk_from_cp, sk_from_cp, mod);
    }
    for(auto& pair: cts)
    {
        pair.first.ntt(ctx);
        pair.second.ntt(ctx);
        GentryPoly::mont_encode(pair.first,  pair.first );
        GentryPoly::mont_encode(pair.second, pair.second);
    }
    std::vector<std::pair<GentryPoly, GentryPoly>> cts_cuda;
    for(auto& pair: cts)
    {
        cts_cuda.push_back({pair.first.to_cuda(), pair.second.to_cuda()});
    }
    return KeySwitchKeyGP(std::move(cts), std::move(cts_cuda), qo);
}
std::pair<GentryPoly, GentryPoly> KeySwitchKeyGP::_key_switch_big_1_cpu(const GentryPoly &a ,const GentryPolyCtx& ctx) const
{
    assert(a.is_cpu());
    std::vector<std::vector<uint64_t>> split = GentryPoly(a).split_by_moduli();

    GentryPoly ra = GentryPoly::zeros_like(cts_[0].first, GPDevice::CPU);
    GentryPoly rb = GentryPoly::zeros_like(cts_[0].first, GPDevice::CPU);
    GentryPoly piece = GentryPoly::zeros_like(cts_[0].first, GPDevice::CPU);

    for(int i=0; i<split.size(); i++)
    {
        auto [cta, ctb] = cts_[i];
        piece.set_from_vec64(split[i]);
        piece.ntt(ctx);
        GentryPoly::mul_mont(cta, cta, piece);
        GentryPoly::mul_mont(ctb, ctb, piece);
        GentryPoly::add(ra, ra, cta);
        GentryPoly::add(rb, rb, ctb);
    }
    ra.intt(ctx);
    rb.intt(ctx);
    ra.moduli_reduce(qo_);
    rb.moduli_reduce(qo_);
    return std::make_pair(std::move(ra), std::move(rb));
}
std::pair<GentryPoly, GentryPoly> KeySwitchKeyGP::_key_switch_big_1_cuda(const GentryPoly &a ,const GentryPolyCtx& ctx) const
{
    assert(a.is_cuda());
    std::vector<CudaBuffer> split = GentryPoly(a).split_by_moduli_cuda();
    GentryPoly ra = GentryPoly::zeros_like(cts_[0].first, GPDevice::CUDA);
    GentryPoly rb = GentryPoly::zeros_like(cts_[0].first, GPDevice::CUDA);
    GentryPoly piece = GentryPoly::zeros_like(cts_[0].first, GPDevice::CUDA);
    
    for(int i=0; i<split.size(); i++)
    {
        auto [cta, ctb] = cts_cuda_[i];
        piece.set_from_cuda_buffer(split[i]);
        piece.ntt(ctx);
        GentryPoly::mul_mont(cta, cta, piece);
        GentryPoly::mul_mont(ctb, ctb, piece);
        GentryPoly::add(ra, ra, cta);
        GentryPoly::add(rb, rb, ctb);
    }
    ra.intt(ctx);
    rb.intt(ctx);
    GentryPoly rac = ra.to_cpu(), rbc = rb.to_cpu();
    rac.moduli_reduce(qo_);
    rbc.moduli_reduce(qo_);
    return std::make_pair(rac.to_cuda(), rbc.to_cuda());
}
std::pair<GentryPoly, GentryPoly> KeySwitchKeyGP::key_switch_big_1(const GentryPoly &a ,const GentryPolyCtx& ctx) const
{
    if(a.is_cpu())return _key_switch_big_1_cpu(a, ctx);
    else return _key_switch_big_1_cuda(a, ctx);
}
std::pair<GentryPoly, GentryPoly> KeySwitchKeyGP::key_switch_big_2(const GentryPoly& a, const GentryPoly& b ,const GentryPolyCtx& ctx) const
{
    auto [ra, rb] = key_switch_big_1(a, ctx);
    GentryPoly::add(rb, rb, b);
    return std::make_pair(std::move(ra), std::move(rb));
}


