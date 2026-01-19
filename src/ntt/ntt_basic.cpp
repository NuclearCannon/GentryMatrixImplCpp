#include <vector>
#include <cassert>
#include <unordered_map>
#include <cmath>
#include "ntt.hpp"


// Cache for bit-reversal tables: key = n (must be power of two)
std::unordered_map<size_t, std::vector<size_t>> bitrev_cache;

const std::vector<size_t>& get_bit_reverse_table(size_t n) {
    auto it = bitrev_cache.find(n);
    if (it != bitrev_cache.end()) {
        return it->second;
    }

    // Compute bits = log2(n)
    size_t bits = 0;
    size_t temp = n;
    while (temp > 1) {
        bits++;
        temp >>= 1;
    }
    assert((1ULL << bits) == n);

    std::vector<size_t> rev(n, 0);
    for (size_t i = 1; i < n; ++i) {
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1ULL) << (bits - 1));
    }

    bitrev_cache[n] = rev;
    return bitrev_cache[n];
}


