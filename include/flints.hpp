#pragma once
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod_vec.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mod_poly.h>
#include <string>
#include <vector>
#include <cassert>

// 各类转化函数

std::string fmpz_to_string(fmpz_t x);
void string_to_fmpz(const std::string s, fmpz_t x);

class fmpz_vector
{
private:
    fmpz* data_;
    int len_;
public:

    fmpz_vector(int len);
    fmpz_vector(const std::vector<std::string>&);

    ~fmpz_vector();

    fmpz_vector(const fmpz_vector&); // 拷贝构造

    fmpz_vector(fmpz_vector&&); // 移动构造

    void print() const;

    inline fmpz* raw(){
        return data_;
    }
    inline const fmpz* raw() const {
        return data_;
    }

    inline fmpz* at(int idx) {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }

    inline fmpz* operator[](int idx) {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }

    inline const fmpz* at(int idx) const {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }

    inline const fmpz* operator[](int idx) const {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }


    std::vector<std::string> export_to_vec_str() const;

};