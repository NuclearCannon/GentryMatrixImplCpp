下面提供Python版本的代码，这会有利于你理解`cuda_ntt.cu`中的`_butterfly_new`函数


```python
import sympy


def find_prime(lower_bound, n):
    if n < 2:
        raise ValueError("n 必须是大于等于 2 的整数")
    if lower_bound < 2:
        lower_bound = 2

    # 从 lower_bound 开始向上找质数
    q = sympy.nextprime(lower_bound - 1)  # 确保包括 lower_bound 本身如果是质数
    while True:
        if q % n == 1:
            return q
        q = sympy.nextprime(q)

def get_q_and_omega_from_n(n):
    q = find_prime(10000, n)
    assert isinstance(q, int)
    qr = sympy.primitive_root(q)
    assert isinstance(qr, int)
    omega = pow(qr, (q-1)//n, q)
    return q, omega

def bit_reverse(x, n, check=True):
    """将x的二进制位反转，用于DIF-NTT的位反转置换"""
    result = 0
    logn = n.bit_length() - 1
    for i in range(logn):
        if (x >> i) & 1:
            result |= 1 << (logn -1 - i)
    if check:
        assert bit_reverse(result, n, check=False) == x
    return result

def bit_reverse_copy(a):
    n = len(a)
    assert (n & (n-1)) == 0
    return [a[bit_reverse(i, n)] for i in range(n)]

def Log2(N):
    assert N & (N-1) == 0
    i = 0
    while 1<<i != N:
        i+=1
    return i


# ==============================================================

def ntt_standard(a, n, q, omega):
    a = a[:]
    # 从大块到小块
    m = n
    t = 0   # t == log2(n/m)
    omegas = [pow(omega, i, q) for i in range(n//2)]
    while m > 1:
        m_half = m>>1
        for k in range(0, n, m):
            for j in range(0, m_half):
                w = omegas[j<<t]    # j<<t == j*n/m <= m_half*n/m = n/2
                i1 = k+j
                i2 = i1 + m_half
                u = a[i1]
                v = a[i2]
                a[i1] = (u + v) % q
                a[i2] = ((u - v) * w) % q
        m = m_half
        t+=1
    # 加上bit reverse就是一个完整的NTT
    return bit_reverse_copy(a)


def ntt_new_version(a, n, q, omega):
    """
    无需bit reverse过程的NTT
    需要(n/2)*log(n) + log(n) + n - 1次乘法
    """
    a = a[:]
    b = a[:]    # 双缓冲区策略
    logn = Log2(n)
    N = n//2    # 记住：N代表n/2
    omegas = [pow(omega, i, q) for i in range(n//2)]    # 可预计算
    for t in range(logn):
        T = 1<<t
        for j in range(0, N, T):
            w = omegas[j]   # 你也可以迭代生成w，但是即使直接读取也只需要n-1次访存，已经很友好了
            for k in range(0, T):                
                i1_r0 = j + k
                i2_r0 = i1_r0 + N
                i1_r1 = 2*j + k
                i2_r1 = i1_r1 + T
                u = a[i1_r0]
                v = a[i2_r0]
                b[i1_r1] = (u + v) % q
                b[i2_r1] = ((u - v) * w) % q   # 这个乘法会被执行log(n)*(n/2)次
        a, b = b, a
    return a


# ==============================================================


def test_ntt_correct(n, q, omega, func):
    a = [i for i in range(n)]
    b1 = func(a, n, q, omega)
    b2 = ntt_standard(a, n, q, omega)
    assert b1 == b2



if __name__ == "__main__":
    n = 256
    q, omega = get_q_and_omega_from_n(n)
    test_ntt_correct(n, q, omega, ntt_new_version)
    # 没assert fail就是对了

    
```