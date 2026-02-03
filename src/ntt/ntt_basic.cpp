#include <vector>
#include <cassert>
#include "ntt.hpp"
#include <memory>

// Cache for bit-reversal tables: key = n (must be power of two)
std::shared_ptr<std::vector<size_t>> bitrev_cache[64];

const std::vector<size_t>& get_bit_reverse_table(size_t n) {
    int log2n = log2(n);
    return get_bit_reverse_table_by_logn(log2n);
}

const std::vector<size_t>& get_bit_reverse_table_by_logn(size_t log2n) {
    if (bitrev_cache[log2n].get())
    {
        return *bitrev_cache[log2n];
    }
    size_t n = 1<<log2n;
    std::shared_ptr<std::vector<size_t>> rev = std::make_shared<std::vector<size_t>>(n, 0);
    for (size_t i = 1; i < n; ++i) {
        (*rev)[i] = ((*rev)[i >> 1] >> 1) | ((i & 1ULL) << (log2n - 1));
    }

    bitrev_cache[log2n] = rev;
    return *rev;
}

