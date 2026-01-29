#pragma once
#include "u64_array.hpp"

std::pair<CRTArray, CRTArray> encrypt64(const CRTArray& message, const CRTArray& sk);

// (m, s, a, b)分别为(系数, NTT, NTT, NTT)形式
// 这种加密适合生成KSK
std::pair<CRTArray, CRTArray> encrypt64_CNNN(const CRTArray& message, const CRTArray& sk);
// 调试用加密函数

// 无噪声
std::pair<CRTArray, CRTArray> encrypt64_no_e(const CRTArray& message, const CRTArray& sk);

// return (0, message)
std::pair<CRTArray, CRTArray> encrypt64_no_ea(const CRTArray& message, const CRTArray& sk);

CRTArray decrypt64(const CRTArray& ct_a, const CRTArray& ct_b, const CRTArray& sk);
