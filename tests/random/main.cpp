#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"
#include "ziq_array.hpp"
#include <iostream>





int main()
{
    fmpz_scalar q("31525889");
    fmpz_vector v0 = fmpz_vector::zeros(20);
    fmpz_vector v1 = fmpz_vector::uniform(20, q.raw());
    fmpz_vector v2 = fmpz_vector::dg(20);
    v0.print();
    v1.print();
    v2.print();
    return 0;
}