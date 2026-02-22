void test_ntt();
void test_ntt_cuda();
void test_encrypt();
void test_encrypt_cuda();
void test_ks();
void test_ks_cuda();
int main()
{
    // test_ntt();
    // test_ntt_cuda();
    // test_encrypt();
    // test_encrypt_cuda();
    test_ks();
    test_ks_cuda();
    return 0;
}