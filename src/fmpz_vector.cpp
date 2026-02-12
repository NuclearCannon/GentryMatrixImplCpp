
#include "flints.hpp"
#include "random.hpp"
#include "uint64.hpp"


std::string fmpz_to_string(fmpz_t x)
{
    char* str = fmpz_get_str(nullptr, 10, x); // 获取十进制字符串
    std::string result(str);
    flint_free(str); // 释放内存
    return result;
}


void string_to_fmpz(const std::string& s, fmpz_t x)
{
    int error =  fmpz_set_str(x, s.c_str(), 10);
    if (error)
    {
        throw std::invalid_argument("Invalid integer string: " + s);
    }
}



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

fmpz_vector::fmpz_vector(const fmpz_vector& other)
{
    int len = other.len_;
    assert(len > 0);
    fmpz *vec = _fmpz_vec_init(len);
    this->data_ = vec;
    this->len_ = len;
    _fmpz_vec_set(vec, other.data_, len);
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

fmpz_vector::fmpz_vector(const std::vector<u64>& src)
{
    int len = src.size();
    fmpz *vec = _fmpz_vec_init(len);
    this->data_ = vec;
    this->len_ = len;
    for(int i=0;i<len;i++)
    {
        fmpz_set_ui(data_+i, src[i]);
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

fmpz_vector fmpz_vector::zeros(int len)
{
    fmpz_vector result(len);
    for(int i=0;i<len;i++)fmpz_zero(result[i]);
    return result;
}
fmpz_vector fmpz_vector::uniform(int len, const fmpz_t q)
{
    // TODO: 可能不够安全！
    fmpz_vector result(len);
    flint_rand_t state;
    flint_randinit(state);  // 初始化随机数生成器
    for (slong i = 0; i < len; i++) {
        fmpz_randm(result[i], state, q);  // vec[i] = random in [0, q)
    }
    flint_randclear(state);  // 清理随机状态
    return result;
}
fmpz_vector fmpz_vector::dg(int len)
{
    // TODO: 这里使用的DG可能不够安全
    fmpz_vector result(len);
    for(int i=0;i<len;i++)fmpz_set_si(result[i], ramdom_generators::dg(10));
    // 
    return result;
}

fmpz_vector fmpz_vector::mod_centered(const fmpz_t q) const
{
    fmpz_vector result(len_);
    fmpz_t q_half, tmp;
    fmpz_init(q_half);
    fmpz_init(tmp);

    // q_half = q / 2
    fmpz_fdiv_q_2exp(q_half, q, 1);

    for (int i = 0; i < len_; ++i) {
        // tmp = data_[i] mod q
        fmpz_mod(tmp, data_ + i, q);

        // if tmp > q/2, tmp -= q
        if (fmpz_cmp(tmp, q_half) > 0) {
            fmpz_sub(tmp, tmp, q);
        }
        fmpz_set(result[i], tmp);
    }

    fmpz_clear(q_half);
    fmpz_clear(tmp);
    return result;
}


// 返回逐位绝对值的最大值，若超出int范围返回-1
long fmpz_vector::max_abs() const
{
    fmpz_t tmp;
    fmpz_init(tmp);
    long max_val = 0;
    for (int i = 0; i < len_; ++i) {
        fmpz_abs(tmp, data_ + i);
        if (fmpz_cmp_si(tmp, max_val) > 0) {
            if (!fmpz_fits_si(tmp)) {
                fmpz_clear(tmp);
                return -1;
            }
            max_val = fmpz_get_si(tmp);
        }
    }
    fmpz_clear(tmp);
    return max_val;
}


std::vector<u64> fmpz_vector::to_uint64() const
{
    std::vector<u64> result(len_);
    static_assert(std::is_same_v<mp_limb_t, u64>);
    for (int i = 0; i < len_; ++i) {
        if (!fmpz_abs_fits_ui(data_+i))
        {
            throw std::invalid_argument("fmpz_vector::to_uint64() 绝对值太大\n");
        }
        if (fmpz_sgn(data_+i) < 0)
        {
            throw std::invalid_argument("fmpz_vector::to_uint64() 存在负数\n");
        }
        result[i] = fmpz_get_ui(data_+i);
    }
    return result;
}