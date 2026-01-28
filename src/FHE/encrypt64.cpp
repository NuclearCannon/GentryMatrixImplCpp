#include "u64_array.hpp"

std::pair<CRTArray, CRTArray> encrypt64(const CRTArray& message, const CRTArray& sk)
{
    auto cc = message.get_cc();
    auto cc_sk = sk.get_cc();
    assert(cc.get() == cc_sk.get());
    CRTArray a = CRTArray::uniform(cc);
    CRTArray e = CRTArray::dg(cc);
    CRTArray as = a.mul_poly(sk);
    CRTArray b = message.add(e).sub(as);
    return std::make_pair(
        std::move(a),
        std::move(b)
    );
}

// 调试用加密函数

// 无噪声
std::pair<CRTArray, CRTArray> encrypt64_no_e(const CRTArray& message, const CRTArray& sk)
{
    auto cc = message.get_cc();
    auto cc_sk = sk.get_cc();
    assert(cc.get() == cc_sk.get());
    CRTArray a = CRTArray::uniform(cc);
    CRTArray as = a.mul_poly(sk);
    CRTArray b = message.sub(as);
    return std::make_pair(
        std::move(a),
        std::move(b)
    );
}

// return (0, message)
std::pair<CRTArray, CRTArray> encrypt64_no_ea(const CRTArray& message, const CRTArray& sk)
{
    auto cc = message.get_cc();
    auto cc_sk = sk.get_cc();
    assert(cc.get() == cc_sk.get());
    return std::make_pair(
        CRTArray::zeros(cc),
        message
    );
}

CRTArray decrypt64(const CRTArray& ct_a, const CRTArray& ct_b, const CRTArray& sk)
{
    // m=as+b
    return ct_a.mul_poly(sk).add(ct_b);
}
