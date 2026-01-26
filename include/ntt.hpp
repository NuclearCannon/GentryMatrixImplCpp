#pragma once
#include <vector>
#include <unordered_map>
#include "flints.hpp"


const std::vector<size_t>& get_bit_reverse_table(size_t n);
int log2(int x);

void ntt_standard_flint(const fmpz* a, fmpz* dst, const fmpz_t root, size_t n, const fmpz_mod_ctx_t ctx);

void ntt_standard_flint(const fmpz_vector& a, fmpz_vector& dst, const fmpz_t root, size_t n, const fmpz_mod_ctx_t ctx);


void ntt_standard_flint_with_roots(
    const fmpz* a, 
    fmpz* dst, 
    const fmpz* roots,   // root的至少[0,n/2)次方 
    size_t n, 
    const fmpz_mod_ctx_t ctx
);


void ntt_standard_flint_with_roots(
    const fmpz_vector& a, 
    fmpz_vector& dst, 
    const fmpz_vector& roots,   // root的至少[0,n/2)次方 
    size_t n, 
    const fmpz_mod_ctx_t ctx
);

void ntt_standard_64(
    const u_int64_t* a, 
    u_int64_t* dst, 
    u_int64_t root,
    size_t n, 
    const u_int64_t mod
);

void ntt_standard_64_with_roots(
    const u_int64_t* a, 
    u_int64_t* dst, 
    const u_int64_t* roots,   // 需要提供root的至少[0,n/2)次方 
    size_t n, 
    const u_int64_t mod
);