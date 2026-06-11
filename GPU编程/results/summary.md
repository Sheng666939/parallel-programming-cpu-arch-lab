# CUDA GPU NTT 实验结果汇总

## 1. 实验环境摘要

- CPU/GPU 平台：Windows 本地笔记本环境。
- GPU：NVIDIA GeForce RTX 4060 Laptop GPU。
- CUDA runtime version：11070。
- CUDA driver version：12060。
- Compute capability：8.9。
- CUDA block size：256。
- 编译器：本机 `nvcc` 为 CUDA 11.7；由于 CUDA 11.7 不支持 `sm_89`，`build.bat` 自动回退到通用 CUDA 编译命令。
- 性能模式状态：benchmark 前后 `nvidia-smi` 均显示 GPU 处于 P0，显存空闲，结果记录在 `results/gpu_status_before.txt` 和 `results/gpu_status_after.txt`。
- 正式 benchmark 命令：`bin\gpu_ntt.exe --min-log 10 --max-log 22 --repeat 20 --warmup 5`。

## 2. 正确性测试结论

已运行：

```powershell
bin\gpu_ntt.exe --verify-only
```

测试覆盖 `poly_n = 1, 2, 3, 7, 16, 31, 127, 513, 1000, 4097`，并覆盖全 0、全 1、递增序列、接近 MOD 的序列、固定随机种子和稀疏多项式。小规模同时对比 O(n^2) 朴素卷积、CPU NTT、GPU naive、GPU Barrett 和 GPU Montgomery；大规模对比 CPU NTT 与三个 GPU 版本。

结论：`ALL VERIFY TESTS PASSED`。

## 3. 各规模时间表

表中 CPU 为串行 CPU NTT baseline；GPU 三列为 end-to-end mean，总时间包含 GPU 函数内的分配、H2D、kernel、D2H 和释放。

| transform_n | poly_n | CPU ms | GPU naive ms | GPU Barrett ms | GPU Montgomery ms |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1024 | 512 | 0.076100 | 0.703540 | 0.771135 | 0.908355 |
| 4096 | 2048 | 0.365800 | 0.786345 | 0.888090 | 0.709195 |
| 16384 | 8192 | 1.379300 | 0.892255 | 0.828745 | 1.216215 |
| 65536 | 32768 | 6.558800 | 1.059205 | 1.015980 | 1.068925 |
| 262144 | 131072 | 30.990300 | 2.520780 | 2.745375 | 2.401610 |
| 1048576 | 524288 | 117.666400 | 8.086115 | 8.209665 | 6.422720 |
| 4194304 | 2097152 | 561.601700 | 28.393265 | 31.660500 | 24.196315 |

## 4. GPU 相对 CPU 加速比

| transform_n | naive | Barrett | Montgomery |
| ---: | ---: | ---: | ---: |
| 1024 | 0.108167 | 0.098686 | 0.083778 |
| 4096 | 0.465190 | 0.411895 | 0.515796 |
| 16384 | 1.545859 | 1.664324 | 1.134092 |
| 65536 | 6.192191 | 6.455639 | 6.135884 |
| 262144 | 12.293933 | 11.288185 | 12.903969 |
| 1048576 | 14.551660 | 14.332668 | 18.320338 |
| 4194304 | 19.779398 | 17.738245 | 23.210216 |

小规模下 GPU end-to-end 时间受 kernel launch、内存分配和拷贝影响明显，因此相对 CPU 不一定加速。从 `2^14` 开始，GPU 版本开始体现并行优势；到 `2^22` 时，Montgomery 版本达到约 23.21 倍 end-to-end 加速。

## 5. Barrett 相对 naive 的提升比例

这里使用 kernel-only 时间计算 `naive_kernel / barrett_kernel`，大于 1 表示 Barrett kernel 更快。

| transform_n | Barrett / naive |
| ---: | ---: |
| 1024 | 0.961131 |
| 4096 | 1.039946 |
| 16384 | 1.025948 |
| 65536 | 1.006727 |
| 262144 | 0.797003 |
| 1048576 | 0.998228 |
| 4194304 | 0.748392 |

Barrett 在中小规模部分测试略优于 naive，但在较大规模中未稳定领先。原因可能是当前 kernel 中 twiddle 仍在线程内计算，整体开销不仅由一次模乘约简决定；同时不同 reduction 路径的寄存器压力和指令调度也会影响最终 kernel 时间。

## 6. Montgomery 相对 naive 的提升比例

这里使用 kernel-only 时间计算 `naive_kernel / montgomery_kernel`，大于 1 表示 Montgomery kernel 更快。

| transform_n | Montgomery / naive |
| ---: | ---: |
| 1024 | 0.695965 |
| 4096 | 1.111208 |
| 16384 | 0.906329 |
| 65536 | 0.973351 |
| 262144 | 1.226269 |
| 1048576 | 1.701096 |
| 4194304 | 1.326664 |

Montgomery 在大规模下表现更好，尤其 `2^20` 和 `2^22`，说明当 butterfly 数量足够多时，Montgomery reduction 的优势能够摊薄域转换成本，并减少普通取模路径的开销。

## 7. kernel-only 与 end-to-end 差异分析

kernel-only 只统计 GPU 计算 kernel，end-to-end 还包含显存分配、H2D/D2H 拷贝、同步和释放。小规模时，kernel 时间本身很短，固定开销占比高，因此 GPU end-to-end 时间可能慢于 CPU。随着规模增大，butterfly 数量快速增加，kernel 并行吞吐成为主导，GPU 相对 CPU 的加速比明显提升。

以 `transform_n = 2^22` 为例，naive 的 kernel-only mean 为 12.228250 ms，而 end-to-end mean 为 28.393265 ms；Montgomery 的 kernel-only mean 为 9.217294 ms，而 end-to-end mean 为 24.196315 ms。可以看到数据传输、分配和同步仍是端到端性能的重要组成部分。

## 8. 报告结论草稿

本实验实现了基于 CUDA C++ 的 NTT 多项式乘法 GPU 并行加速。CPU baseline 采用迭代 Cooley-Tukey NTT，GPU 版本采用每个线程负责一个 butterfly、每个 stage 独立 kernel 的并行策略，并实现了 naive、Barrett 和 Montgomery 三种模乘路径。正确性测试覆盖多种规模和数据分布，所有 GPU 结果均与 CPU NTT reference 一致。性能实验表明，小规模输入受 kernel launch、内存分配和 H2D/D2H 拷贝影响较大，GPU end-to-end 性能不一定优于 CPU；当 transform size 增大后，GPU 并行 butterfly 的吞吐优势逐渐显现。在 `2^22` 规模下，Montgomery 版本取得约 23.21 倍 end-to-end 加速，说明对大量模乘场景，使用 Montgomery reduction 具有较明显的优化价值。后续可继续通过 twiddle 预处理、shared memory stage 合并、batch NTT 和 CUDA Graph 降低重复计算与 kernel launch 开销。
