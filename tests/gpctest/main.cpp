void test_ntt();
void test_ntt_cuda();
void test_encrypt();
void test_encrypt_cuda();
void test_ks();
int main()
{
    // test_ntt();
    // test_ntt_cuda();
    test_encrypt();
    test_encrypt_cuda();
    // test_ks();
    return 0;
}