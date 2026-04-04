
#include "flints.hpp"
#include "random.hpp"
#include <cstdint>



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

fmpz_vector::fmpz_vector(const std::vector<uint64_t>& src)
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
    for(int i=0;i<len;i++)fmpz_set_si(result[i], random_generators::dg(10));
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
int64_t fmpz_vector::max_abs() const
{
    fmpz_t tmp;
    fmpz_init(tmp);
    int64_t max_val = 0;
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


std::vector<uint64_t> fmpz_vector::to_uint64() const
{
    std::vector<uint64_t> result(len_);
    static_assert(std::is_same_v<mp_limb_t, uint64_t>);
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



// CRT分解
std::vector<std::vector<uint64_t>> crt(const fmpz_vector& src, const std::vector<uint64_t>& mods)
{
    std::vector<std::vector<uint64_t>> result;
    int len = src.len();
    fmpz_vector buf(len);
    for(uint64_t mod : mods)
    {
        fmpz_scalar mod_mpz = fmpz_scalar::from_ui(mod);
        _fmpz_vec_scalar_mod_fmpz(buf.raw(), src.raw(), len, mod_mpz.raw());
        result.push_back(buf.to_uint64());
    }
    return result;
}

// CRT合并
void icrt(fmpz_vector& dst, const std::vector<std::vector<uint64_t>>& src, const std::vector<uint64_t>& mods)
{
    int mod_len = mods.size();
    assert(mod_len >= 1);
    int data_len = src[0].size();

    // 初始化结果向量 x 为第一个模下的值
    fmpz_vector x(data_len);
    for (int j = 0; j < data_len; ++j) {
        fmpz_set_ui(x[j], src[0][j]);  // x_j ≡ src[0][j] (mod mods[0])
    }

    // 当前累积模数 M = mods[0]
    fmpz_scalar M;
    fmpz_set_ui(M.raw(), mods[0]);

    // 临时变量用于 fmpz_crt
    fmpz_t r1, r2, m1, m2, res;
    fmpz_init(r1); fmpz_init(r2);
    fmpz_init(m1); fmpz_init(m2);
    fmpz_init(res);

    // 逐个合并后续模数
    for (int i = 1; i < mod_len; ++i) {
        uint64_t mi = mods[i];
        fmpz_set_ui(m2, mi);          // m2 = mods[i]
        fmpz_set(m1, M.raw());        // m1 = 当前累积模数 M

        for (int j = 0; j < data_len; ++j) {
            fmpz_set(r1, x[j]);               // r1 = 当前解 x_j
            fmpz_set_ui(r2, src[i][j]);       // r2 = 新余数 src[i][j]

            // 调用 FLINT 的 CRT：res ≡ r1 (mod m1), res ≡ r2 (mod m2)
            fmpz_CRT(res, r1, m1, r2, m2, 0);  // 0 => 非负最小解

            fmpz_set(x[j], res);              // 更新 x[j]
        }

        // 更新累积模数 M = M * mods[i]
        fmpz_mul_ui(M.raw(), M.raw(), mi);
    }

    // 写入输出
    assert(dst.len() == data_len);
    _fmpz_vec_set(dst.raw(), x.raw(), data_len);

    // 清理临时变量
    fmpz_clear(r1); fmpz_clear(r2);
    fmpz_clear(m1); fmpz_clear(m2);
    fmpz_clear(res);
}
