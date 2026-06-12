# GPU NTT 多项式乘法实验记录草稿

本实验环境中 CUDA runtime version 为 11070，即 CUDA 11.7；CUDA driver version 为 12060，即 CUDA 12.6。

## 1. NTT 并行化思路

NTT 的每一层 stage 由若干个互不重叠的 butterfly 组成。在同一 stage 内，不同 butterfly 访问的系数位置互不冲突，因此可以让每个 CUDA thread 负责一个 butterfly。

当前 GPU kernel 的线程映射为：

```text
tid = blockIdx.x * blockDim.x + threadIdx.x
group = tid / mid
offset = tid % mid
pos1 = group * 2 * mid + offset
pos2 = pos1 + mid
```

每个线程读取 `a[pos1]` 和 `a[pos2]`，计算 `v = a[pos2] * w[offset] mod MOD`，再写回加法和减法结果。

## 2. Stage 之间的依赖

同一 stage 内的 butterfly 可以并行执行，但 stage 之间存在严格的数据依赖：后一层必须使用前一层已经完成的结果。CUDA block 之间没有直接的全局同步机制，因此当前实现为每个 stage 单独启动一个 kernel，利用 kernel launch 之间的天然全局同步保证正确性。

## 3. GPU 版本主要瓶颈

- 模乘开销：朴素版本使用 64 位乘法后 `% MOD`，设备端整数除法代价较高。
- 访存 stride：随着 stage 增大，butterfly 的两个输入位置跨度增大，访存局部性会变化。
- kernel launch overhead：每个 stage 一个 kernel，小规模 NTT 中 launch 开销占比明显。
- H2D/D2H 拷贝：端到端时间包含数据传输，小规模时拷贝和同步可能掩盖 kernel 加速效果。
- twiddle 计算：当前第一阶段在线程内用快速幂计算 `w[offset]`，实现简单但会增加算术开销。

## 4. Barrett 和 Montgomery 优化原因

Barrett reduction 通过预计算 `floor(2^64 / MOD)`，使用高 64 位乘法近似商，避免直接执行 `% MOD`，因此可能减少模乘中的整数除法开销。

Montgomery reduction 将数转换到 Montgomery 域中，模乘时使用低位乘法和移位完成约简。对于需要大量模乘的 NTT，若域转换成本能被足够多的 butterfly 摊薄，Montgomery 版本可能优于朴素取模版本。

## 5. 后续优化方向

- twiddle 预处理：将每个 stage 的旋转因子提前写入 GPU 内存，避免每个线程重复快速幂。
- shared memory 小 stage 合并：对较小 stage 使用 shared memory，并在一个 block 内合并多层计算，减少 global memory 访问和 kernel launch 次数。
- batch NTT：同时处理多个多项式，提升 GPU 占用率并摊薄 launch overhead。
- pinned memory 与异步拷贝：优化 H2D/D2H 拷贝时间。
- 使用 CUDA Graph：固定 benchmark 流程时可减少重复 kernel launch 管理开销。
