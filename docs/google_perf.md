
```
# 安装（Ubuntu/Debian）
sudo apt install libgoogle-perftools-dev google-perftools

# 编译时链接 tcmalloc 和 profiler
g++ -o bench bench.cpp -lprofiler -ltcmalloc

# 运行并生成 profile
./bench

# 查看结果
pprof --text ./bench bench.prof      # 文本
pprof --pdf ./bench bench.prof > out.pdf  # PDF 调用图
pprof --web ./bench bench.prof # 网页
pprof -kcachegrind ./bench bench.prof # 使用kcachegrind打开

```