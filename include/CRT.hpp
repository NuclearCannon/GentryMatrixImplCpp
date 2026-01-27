#pragma once
#include "flints.hpp"
#include <vector>
#include "uint64.hpp"

// CRT分解
vv64 crt(const fmpz_vector& src, const vec64& mods);

// CRT合并
void icrt(fmpz_vector& dst, const vv64& src, const vec64& mods);
