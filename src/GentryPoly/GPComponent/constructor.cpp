#include "GentryPoly.hpp"

GPComponent::GPComponent(size_t n, size_t p, uint64_t q):
    n_(n),
    p_(p),
    q_(q),
    mm_(q),
    data_(2*n*n*(p-1))
{

}

// 私有裸构造
GPComponent::GPComponent(size_t n, size_t p, uint64_t q, GPComponent::tag_no_data):
    n_(n),
    p_(p),
    q_(q),
    mm_(q)
{
}

GPComponent GPComponent::from_data(size_t n, size_t p, uint64_t q, std::vector<uint64_t> data)
{
    // data形参是值类型而不是引用，因为我们需要一个完整对象所有权
    GPComponent result(n, p, q, tag_no_data {});
    // 检查data长度
    assert(data.size() == 2*n*n*(p-1));
    // TODO: 要不要检查data取值？
    result.data_.swap(data);
    return result;
}


GPComponent::~GPComponent() = default;
GPComponent::GPComponent(const GPComponent&) = default;
GPComponent::GPComponent(GPComponent&&) = default;