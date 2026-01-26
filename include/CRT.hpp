#pragma once
#include "flints.hpp"
#include <vector>
#include <ctype.h>

// CRT分解
std::vector<std::vector<u_int64_t>> crt(const fmpz_vector& src, const std::vector<u_int64_t>& mods);

// CRT合并
void icrt(fmpz_vector& dst, const std::vector<std::vector<u_int64_t>>& src, const std::vector<u_int64_t>& mods);
