#include <cstring>
#include <cstdio>

void bench_ks();
void bench_ks_cuda();
void bench_ntt_cuda();
void bench_circledast();
void bench_circledast_cuda();

int main(int argc, char *argv[])
{
    for(int i=1; i<argc; i++)   // 从1开始以忽略ELF文件名
    {
        if(strcmp(argv[i], "--ks")==0)bench_ks();
        else if(strcmp(argv[i], "--ksc")==0)bench_ks_cuda();
        else if(strcmp(argv[i], "--nttc")==0)bench_ntt_cuda();
        else if(strcmp(argv[i], "--cd")==0)bench_circledast();
        else if(strcmp(argv[i], "--cdc")==0)bench_circledast_cuda();
        else {
            printf("未被识别的选项：%s\n", argv[i]);
        }
    }
    return 0;
}