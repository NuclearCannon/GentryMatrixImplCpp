#include "ziq_array.hpp"

// TODO: 或许调用者可以指定message和sk的格式，而不仅限于系数

std::pair<ZiqArray, ZiqArray> encrypt(const ZiqArray& message, const ZiqArray& sk);


ZiqArray decrypt(const ZiqArray& ct_a, const ZiqArray& ct_b, const ZiqArray& sk);