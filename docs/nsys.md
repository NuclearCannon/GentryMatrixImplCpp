nsys的一种比较通用的用法是：

```
nsys profile -t cuda,osrt,nvtx -o my_app_profile --force-overwrite true --stats true ./my_cuda_app arg1 arg2
```

对于当前版本来说，如果你想要测GPU KS效率
```
nsys profile --capture-range=cudaProfilerApi -t cuda,osrt,nvtx -o my_app_profile --force-overwrite true --stats true ./bench --ksc
```

`--capture-range=cudaProfilerApi`是为了指定我们的监控区间

你可以如此阅读测试结果：

使用`nsys-ui`打开`.nsys-rep`文件

时间轴 (Timeline)：你会看到多行轨道。

找到 CUDA API 行：你会看到密密麻麻的黄色/红色条，那就是 cudaMalloc 和 cudaMemcpy。你会发现它们把时间轴填满了，而绿色的 CUDA Kernel 条被挤得很散。

- 找到 CUDA Kernel 行：观察绿色条之间是否有巨大的空白？那些空白就是 GPU 在等 CPU 分配内存或传数据。

- 找到 OS Runtime 行：你会看到 CPU 线程在 poll 或 futex 上阻塞，对应着等待 GPU 的时刻。

- 悬浮提示：鼠标悬停在任何一个条上，会显示详细的耗时、参数（比如 memcpy 的大小）、调用栈。

- 范围缩放：使用鼠标滚轮或顶部缩放条，放大到具体的某一次 malloc + memcpy + kernel 过程，你会清晰地看到串行执行的低效。


