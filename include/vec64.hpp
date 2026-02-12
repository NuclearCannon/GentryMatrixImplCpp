#pragma once
#include <vector>
#include <type_traits>
#include <cstdint>
#include <cassert>

typedef std::vector<uint64_t> vec64;
typedef std::vector<vec64> vv64;

// 返回x的[0,len)次幂组成的向量
void get_powers(vec64& dst, uint64_t x, size_t len, uint64_t mod);

// 返回x的[0,len)次幂组成的向量
vec64 get_powers(uint64_t x, size_t len, uint64_t mod);


// 不一定有用，但是调试的时候或许很有用。
std::size_t vec64hash(const vec64& v);
std::size_t vv64hash(const vv64& v);
