# CUDA GPU NTT 工程代码审计记录

本文档记录第二阶段对当前 CUDA C++ NTT 多项式乘法工程的代码审计结论。审计目标是确认常量、CPU baseline、GPU kernel、计时口径和 benchmark 流程符合实验报告与 GitHub 提交要求。

## 1. `include/mod_arith.cuh`

### 常量一致性

- 设备端模数使用 `MOD_DEVICE = 998244353u`，主机端公共常量在 `include/common.hpp` 中为 `MOD = 998244353u`，二者一致。
- 原根使用 `G_DEVICE = 3u`，主机端公共常量为 `G = 3u`，二者一致。
- Barrett 常量 `BARRETT_MU = 18479187002ull`，等于 `floor(2^64 / 998244353)`，常量正确。
- Montgomery 常量：
  - `MONT_N_PRIME = 998244351u`，等于 `-MOD^{-1} mod 2^32`。
  - `MONT_R = 301989884u`，等于 `2^32 mod MOD`。
  - `MONT_R2 = 932051910u`，等于 `(2^32)^2 mod MOD`。

### 模运算实现结论

- `add_mod` 和 `sub_mod` 输入输出均保持在 `[0, MOD)`。
- `mul_mod_naive` 使用 64 位乘法后 `% MOD`，作为 baseline 正确。
- `barrett_reduce` 使用 `__umul64hi(x, BARRETT_MU)` 估计商，并进行两次修正，满足本工程中 `x < MOD^2` 的乘法约简需求。
- `mont_reduce` 按 32 位 Montgomery reduction 实现，`mont_mul`、`to_mont_value`、`from_mont_value` 组合正确。
- `bin\gpu_ntt.exe --verify-only` 已覆盖 naive、Barrett、Montgomery 三种设备端模乘路径，所有测试通过，说明三种模乘在当前 NTT 使用路径中结果一致。

## 2. `src/ntt_cpu.cpp`

### bit reversal

- `ntt_cpu` 使用迭代式 bit-reversal 交换逻辑，循环变量 `i, j, bit` 的更新方式为常见 Cooley-Tukey NTT 写法。
- 仅在 `i < j` 时交换，避免重复交换。

### 正变换、逆变换与 `inv_n`

- 每层 `len` 从 2 倍增到 `n`，`mid = len / 2`。
- 正变换 stage root 使用 `G^((MOD - 1) / len)`。
- 逆变换 stage root 使用正向 stage root 的模逆。
- 逆变换结束后统一乘 `inv_n = n^(MOD - 2) mod MOD`，符合 NTT 逆变换要求。
- 所有乘法路径使用 `uint64_t` 中间值，避免 32 位溢出。

### 卷积输出长度

- `convolution_cpu` 的输出长度为 `a.size() + b.size() - 1`。
- 内部 transform 长度使用不小于结果长度的最小 2 的幂。
- 小规模测试中 CPU NTT 与 O(n^2) 朴素卷积结果一致。

## 3. `src/ntt_cuda.cu`

### GPU 核心计算位置

- GPU 版本在设备端完成 bit reversal、每层 NTT butterfly、pointwise multiply、逆 NTT 后乘 `inv_n`。
- Host 端只负责数据准备、stage root 计算、kernel launch、计时与结果拷贝，不承担 NTT butterfly 主体计算。

### stage 同步语义

- 每个 stage 单独启动一次 `ntt_stage_kernel`。
- CUDA kernel launch 顺序提供 stage 之间的全局同步语义，保证后一层读取的是前一层完成后的数据。
- 同一 stage 内每个线程负责一个 butterfly，线程映射为 `tid -> group/offset -> pos1/pos2`，同一 stage 内写入位置互不重叠。

### CUDA 错误检查

- `cudaEventCreate`、`cudaMalloc`、`cudaMemcpy`、`cudaEventRecord`、`cudaEventSynchronize`、`cudaEventElapsedTime` 均通过 `CUDA_CHECK_LOCAL` 检查。
- 每次 kernel launch 后调用 `CUDA_KERNEL_CHECK_LOCAL()`，可及时发现 launch 参数和设备端执行错误。
- 若 CUDA API 失败，函数返回 `correct=false`，benchmark 层继续处理其他算法。

### kernel-only 计时口径

- `kernel_ms` 使用 CUDA event 包围 GPU 计算段。
- kernel-only 包含 bit reversal、NTT stages、pointwise multiply、逆 NTT、Montgomery 域转换等 GPU kernel。
- kernel-only 不包含 CSV 文件 I/O、CPU 正确性比较和 benchmark 汇总输出。
- `total_ms` 使用主机端 `std::chrono`，包含一次 GPU 卷积函数内的内存分配、拷贝、kernel、同步、结果返回和释放。

## 4. `src/benchmark.cpp` / `src/main.cu`

### warmup 与 repeat

- `warmup` 调用 `convolution_cuda` 但不计入正式平均时间。
- `repeat` 循环中记录 H2D、kernel、D2H、total，并输出 mean；`gpu_total_ms_median` 使用 repeat 样本计算中位数。
- CPU baseline 当前每个规模运行一次，因此 `cpu_ms_mean` 与 `cpu_ms_median` 相同。

### CPU baseline 适用范围

- CPU baseline 是同一台机器上的串行 NTT 对照，用于本地 GPU 加速比分析。
- 它不是跨机器绝对性能标准；报告中应说明 CPU/GPU 均来自本机环境。

### 正确性增强

- 已支持 `bin\gpu_ntt.exe --verify-only`。
- 覆盖 `poly_n = 1, 2, 3, 7, 16, 31, 127, 513, 1000, 4097`。
- 覆盖全 0、全 1、递增序列、接近 MOD 的序列、固定随机种子、稀疏多项式。
- 小规模对比 O(n^2)、CPU NTT 和三个 GPU 版本；大规模至少对比 CPU NTT 和三个 GPU 版本。

## 5. 审计结论

当前工程满足第二阶段要求：核心算法在 GPU 上执行，CPU baseline 正确，三种 GPU 模乘路径均通过扩展正确性测试，计时口径区分 kernel-only 与 end-to-end，CSV 字段满足正式 benchmark 记录要求。后续性能优化可以继续围绕 twiddle 预处理、shared memory stage 合并、batch NTT 和 CUDA Graph 展开，但这些优化不属于本阶段提交范围。
