#include "application.hpp"
#include <iostream>

SecretKey::SecretKey(size_t n, size_t p, const std::vector<uint64_t>& mods):
    data_(std::make_unique<GentryPoly>(
        GentryPoly::sk(n, p, mods)
    ))
{

}