#pragma once
#include <vector>
#include <unordered_map>
#include "flints.hpp"


const std::vector<size_t>& get_bit_reverse_table(size_t n);

void ntt_standard_flint(const fmpz* a, fmpz* dst, const fmpz_t root, size_t n, const fmpz_mod_ctx_t ctx);

void ntt_standard_flint(const fmpz_vector& a, fmpz_vector& dst, const fmpz_t root, size_t n, const fmpz_mod_ctx_t ctx);