void test_ntt();
void test_ntt_cuda();
void test_encrypt();
int main()
{
    test_ntt();
    test_ntt_cuda();
    test_encrypt();
    return 0;
}