#pragma once
#include <vector>
#include <unordered_map>
#include "flints.hpp"


extern std::unordered_map<size_t, std::vector<size_t>> bitrev_cache;

const std::vector<size_t>& get_bit_reverse_table(size_t n);

void ntt_standard_flint(const fmpz* a, fmpz* dst, fmpz_t root, size_t n, const fmpz_t q);

void ntt_standard_flint(const fmpz_vector& a, fmpz_vector& dst, fmpz_t root, size_t n, const fmpz_t q);