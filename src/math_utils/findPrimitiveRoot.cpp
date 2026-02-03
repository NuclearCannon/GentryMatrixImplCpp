#include <vector>
#include "uint64.hpp"


// 寻找素数 p 的一个原根
u64 findPrimitiveRoot(u64 p) {
    if (p == 2) return 1;

    u64 phi = p - 1; // φ(p) = p - 1
    
    u64 n = phi;
    std::vector<u64> factors;
    for (u64 i = 2; i * i <= n; ++i) {
        if (n % i == 0) {
            factors.push_back(i);
            while (n % i == 0)
                n /= i;
        }
    }
    if (n > 1)
        factors.push_back(n);

    // 尝试从 2 开始逐个检查是否为原根
    for (u64 g = 2; g < p; ++g) {
        bool ok = true;
        for (u64 q : factors) {
            // 检查 g^(φ/q) ≡ 1 (mod p) 是否成立
            if (mod_pow(g, phi / q, p) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
    return -1; // 理论上不会发生（素数必有原根）
}