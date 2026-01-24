#pragma once
#include "ziq_array.hpp"

// TODO: 或许调用者可以指定message和sk的格式，而不仅限于系数

std::pair<ZiqArray, ZiqArray> encrypt(const ZiqArray& message, const ZiqArray& sk);

// 调试用加密函数

// 无噪声
std::pair<ZiqArray, ZiqArray> encrypt_no_e(const ZiqArray& message, const ZiqArray& sk);

// return (0, message)
std::pair<ZiqArray, ZiqArray> encrypt_no_ea(const ZiqArray& message, const ZiqArray& sk);

ZiqArray decrypt(const ZiqArray& ct_a, const ZiqArray& ct_b, const ZiqArray& sk);

// 加密，(m, sk, cta, ctb)的格式分别为(Coeff, NTT, NTT, NTT)
// 这种加密函数有利于构造KSK, KSK需要NTT形式的密文
std::pair<ZiqArray, ZiqArray> encrypt_CNNN(const ZiqArray& message, const ZiqArray& sk_ntt);