#include "flints.hpp"

fmpz_vector::fmpz_vector(int len)
{
    assert(len > 0);
    fmpz *vec = _fmpz_vec_init(len);
    // 不初始化
    this->data_ = vec;
    this->len_ = len;
}


fmpz_vector::~fmpz_vector()
{
    _fmpz_vec_clear(data_, len_);
    len_ = 0;
}


fmpz_vector::fmpz_vector(fmpz_vector&& other)
{
    data_ = other.data_;
    len_ = other.len_;
    other.data_ = nullptr;
    other.len_ = 0;
}


fmpz_vector::fmpz_vector(const std::vector<std::string>& lst_str)
{
    int len = lst_str.size();
    fmpz *vec = _fmpz_vec_init(len);
    this->data_ = vec;
    this->len_ = len;
    for(int i=0;i<len;i++)
    {
        string_to_fmpz(lst_str[i], vec+i);
    }
}


std::vector<std::string> fmpz_vector::export_to_vec_str() const
{
    std::vector<std::string> result;
    for(int i=0;i<len_;i++)
    {
        result.push_back(
            fmpz_to_string(data_ + i)
        );
    }
    return result;
}

void fmpz_vector::print() const
{
    auto vs = export_to_vec_str();
    printf("[");
    for(auto& s : vs)
    {
        printf("%s, ", s.c_str());
    }
    printf("]\n");
}