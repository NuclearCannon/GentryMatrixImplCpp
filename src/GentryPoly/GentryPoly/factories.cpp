#include "GentryPoly.hpp"
#include "random.hpp"


GentryPoly GentryPoly::zeros(size_t n, size_t p, const std::vector<uint64_t>& moduli)
{
    std::vector<GPComponent> cpu_comps;
    for(uint64_t q : moduli)
    {
        cpu_comps.push_back(GPComponent::zeros(n, p, q));
    }
    return GentryPoly(false, moduli, cpu_comps);
}
GentryPoly GentryPoly::dg(size_t n, size_t p, const std::vector<uint64_t>& moduli)
{
    size_t size = 2*(p-1)*n*n;
    std::vector<int64_t> data(size);
    for(auto& i: data) i = random_generators::dg(5);
    std::vector<GPComponent> cpu_comps;
    for(uint64_t q : moduli)
    {
        cpu_comps.push_back(GPComponent::from_signed_data(n, p, q, data));
    }
    return GentryPoly(false, moduli, cpu_comps);
}
GentryPoly GentryPoly::sk(size_t n, size_t p, const std::vector<uint64_t>& moduli)
{
    size_t size = 2*(p-1)*n*n;
    std::vector<int64_t> data(size, 0);
    for(int i=0;i<size;i+=n)data[i] = random_generators::randint(-1, 1);
    std::vector<GPComponent> cpu_comps;
    for(uint64_t q : moduli)
    {
        cpu_comps.push_back(GPComponent::from_signed_data(n, p, q, data));
    }
    return GentryPoly(false, moduli, cpu_comps);
}
GentryPoly GentryPoly::uniform(size_t n, size_t p, const std::vector<uint64_t>& moduli)
{
    size_t size = 2*(p-1)*n*n;
    std::vector<GPComponent> cpu_comps;
    for(uint64_t q : moduli)
    {
        std::vector<uint64_t> data(size);
        for(auto& i:data)i = random_generators::randu64(0, q-1);
        cpu_comps.push_back(GPComponent::from_data(n, p, q, std::move(data)));
    }
    return GentryPoly(false, moduli, cpu_comps);
}
GentryPoly GentryPoly::randint(size_t n, size_t p, const std::vector<uint64_t>& moduli, int lb, int ub)
{
    size_t size = 2*(p-1)*n*n;
    std::vector<int64_t> data(size);
    for(auto& i: data) i = random_generators::randint(lb, ub);
    std::vector<GPComponent> cpu_comps;
    for(uint64_t q : moduli)
    {
        cpu_comps.push_back(GPComponent::from_signed_data(n, p, q, data));
    }
    return GentryPoly(false, moduli, cpu_comps);
}

GentryPoly GentryPoly::to_cpu() const {
    if(is_cpu())return *this;
    assert(is_cuda());
    GentryPoly res = zeros_like(*this);
    assert(res.is_cpu());
    auto& dst = res.cpu_components();
    auto& src = this->cuda_components();
    int len = dst.size();
    assert(src.size() == len);
    for(int i=0; i<len; i++)
    {
        src[i].to_cpu(dst[i]);
    }
    return res;

}
GentryPoly GentryPoly::to_cuda() const {
    if(is_cuda())return *this;
    assert(is_cpu());
    // 构造Cuda Storage
    CudaStorage sto;
    for(auto& i : cpu_components())
    {
        sto.emplace_back(GPComponentCuda::copy_from_cpu(i));
    }
    return GentryPoly(true, moduli_, sto);
}