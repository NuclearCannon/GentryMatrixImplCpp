#include "application.hpp"
#include <iostream>
#include "random.hpp"

SecretKey::SecretKey(size_t n, size_t p):
    n_(n), p_(p), sk_data_(2*n_*(p_-1))
{
    for(int i=0; i<sk_data_.size(); i++)
    {
        sk_data_[i] = random_generators::randint(-1, 1);
    }
}

GentryPoly SecretKey::as_poly(const std::vector<uint64_t>& mods) const
{
    fmpz_vector tmp(sk_data_.size()*n_);
    for(int i=0; i<sk_data_.size(); i++)
    {
        fmpz_set_si(tmp[i*n_], sk_data_[i]);
    }
    return GentryPoly::from_fmpz(n_, p_, mods, tmp);
}