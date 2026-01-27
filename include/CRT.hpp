#pragma once
#include "flints.hpp"
#include <vector>
#include "uint64.hpp"

// CRT分解
std::vector<std::vector<u64>> crt(const fmpz_vector& src, const std::vector<u64>& mods);

// CRT合并
void icrt(fmpz_vector& dst, const std::vector<std::vector<u64>>& src, const std::vector<u64>& mods);
