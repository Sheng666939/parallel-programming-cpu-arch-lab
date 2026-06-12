# CUDA C++ NTT 多项式乘法 GPU 加速实验

本实验主题是使用 CUDA C++ 实现 NTT 多项式乘法的 GPU 并行加速，模数使用 `998244353`，原根使用 `3`。工程包含 CPU 串行 baseline、GPU naive、GPU Barrett 和 GPU Montgomery 四个版本，并输出 correctness 与性能统计 CSV。

## 实验环境

- 操作系统：Windows
- GPU：本地 NVIDIA RTX 4060
- CUDA runtime version：11070，即 CUDA 11.7
- CUDA driver version：12060，即 CUDA 12.6
- 编译工具：CUDA Toolkit、MSVC、`nvcc`
- 建议终端：`x64 Native Tools Command Prompt for VS`，或已经配置好 MSVC 环境变量的 PowerShell

检查环境：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/check_env.ps1
```

如果 `nvcc` 不存在，需要先安装 CUDA Toolkit，并确认终端中能找到 MSVC 编译器。CUDA 11.7 等较旧版本可能不支持 `sm_89`，`build.bat` 会先尝试 `sm_89`，失败后自动使用通用 CUDA 编译命令。

本机 nvcc 为 CUDA 11.7，无法直接使用 RTX 4060 对应的 sm_89 编译目标，因此 build.bat 自动回退到通用 CUDA 编译命令；本实验结果可以反映当前本地环境下的真实性能，但后续若升级 CUDA 11.8/12.x 并使用 sm_89 重新编译，可能进一步改善性能。

## 编译方式

```bat
build.bat
```

等价核心命令：

```bat
nvcc -O3 -std=c++17 -arch=sm_89 ^
  src/main.cu src/ntt_cpu.cpp src/ntt_cuda.cu src/benchmark.cpp ^
  -Iinclude ^
  -o bin/gpu_ntt.exe
```

如果当前 CUDA Toolkit 不支持 `sm_89`，可以使用：

```bat
nvcc -O3 -std=c++17 src/main.cu src/ntt_cpu.cpp src/ntt_cuda.cu src/benchmark.cpp -Iinclude -o bin/gpu_ntt.exe
```

## 运行方式

默认测试 `2^10, 2^12, ..., 2^22`：

```powershell
powershell -ExecutionPolicy Bypass -File run_bench.ps1
```

也可以直接运行：

```powershell
bin/gpu_ntt.exe --min-log 10 --max-log 22 --repeat 20 --warmup 3
```

小规模快速检查：

```powershell
bin/gpu_ntt.exe --min-log 10 --max-log 10 --repeat 1 --warmup 0
```

扩展正确性验证：

```powershell
bin/gpu_ntt.exe --verify-only
```

正式 benchmark：

```powershell
nvidia-smi --query-gpu=timestamp,name,pstate,power.draw,power.limit,clocks.gr,clocks.mem,temperature.gpu,utilization.gpu,memory.used,memory.total --format=csv > results/gpu_status_before.txt
bin/gpu_ntt.exe --verify-only
bin/gpu_ntt.exe --min-log 10 --max-log 22 --repeat 20 --warmup 5
nvidia-smi --query-gpu=timestamp,name,pstate,power.draw,power.limit,clocks.gr,clocks.mem,temperature.gpu,utilization.gpu,memory.used,memory.total --format=csv > results/gpu_status_after.txt
```

## 代码结构

- `include/common.hpp`：公共常量、benchmark 参数、主机端辅助函数。
- `include/ntt_cpu.hpp`：CPU NTT 与朴素卷积接口。
- `include/ntt_cuda.cuh`：GPU 卷积接口、算法枚举和结果结构。
- `include/mod_arith.cuh`：CUDA 设备端模加、模减、朴素模乘、Barrett、Montgomery 运算。
- `include/timer.hpp`：主机端计时器。
- `src/ntt_cpu.cpp`：CPU 迭代 Cooley-Tukey NTT 实现。
- `src/ntt_cuda.cu`：GPU NTT kernels 和 `convolution_cuda`。
- `src/benchmark.cpp`：数据生成、正确性检查、性能统计和 CSV 写入。
- `src/main.cu`：命令行参数解析。
- `results/bench.csv`：benchmark 输出结果。
- `docs/experiment_notes.md`：实验记录草稿。

## 已实现算法版本

1. CPU baseline：迭代 NTT，多项式卷积使用 64 位乘法防止溢出。
2. GPU naive：每个 CUDA 线程负责一个 butterfly，每个 stage 单独启动 kernel，模乘使用普通 `% MOD`。
3. GPU Barrett：线程映射与 naive 相同，模乘使用 Barrett reduction，减少整数除法开销。
4. GPU Montgomery：输入、twiddle 和中间结果进入 Montgomery 域，使用 Montgomery reduction 进行模乘，最后转换回普通整数域。

## 计时口径

- `gpu_h2d_ms`：Host to Device 拷贝时间。
- `gpu_kernel_ms`：只统计 GPU 计算 kernel，包括 NTT、pointwise multiply、逆 NTT，以及 Montgomery 域转换 kernel。
- `gpu_d2h_ms`：Device to Host 拷贝时间。
- `gpu_total_ms`：从 `convolution_cuda` 开始到结果返回结束的 end-to-end 时间，包含分配、拷贝、kernel、同步和释放。

## CSV 字段说明

`results/bench.csv` 字段包括：

```text
timestamp,gpu_name,cuda_runtime_version,cuda_driver_version,compute_capability,
transform_n,poly_n,algo,repeat,warmup,block_size,
cpu_ms_mean,cpu_ms_median,gpu_h2d_ms_mean,gpu_kernel_ms_mean,
gpu_d2h_ms_mean,gpu_total_ms_mean,gpu_total_ms_median,
speedup_vs_cpu,speedup_vs_gpu_naive,throughput_mbutterfly_per_s,correct
```

吞吐率按照一次多项式乘法约 `3 * (n / 2) * log2(n)` 个 butterfly 计算，并除以 kernel-only 时间。

## 报告分析建议

- 比较 CPU baseline 与 GPU 版本在不同规模下的加速比。
- 比较 GPU naive、Barrett、Montgomery 的 kernel-only 时间。
- 分析 kernel-only 与 end-to-end 时间的差距，观察 H2D/D2H 拷贝对小规模测试的影响。
- 观察规模增大后加速比变化：小规模可能受 kernel launch overhead 主导，大规模更能体现 GPU 并行吞吐。
- 后续可加入 twiddle 预处理、shared memory 小 stage 合并和 batch NTT，对比进一步优化效果。

## 第二阶段结果文件

- `docs/code_audit.md`：代码审计记录，覆盖模运算、CPU NTT、GPU kernel、计时和 benchmark 口径。
- `results/summary.md`：正式 benchmark 汇总，可直接用于实验报告。
- `results/gpu_status_before.txt` / `results/gpu_status_after.txt`：正式 benchmark 前后的 GPU 状态。
