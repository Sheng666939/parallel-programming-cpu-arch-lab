# GPU NTT 实验报告结论草稿

本实验首先完成了 CPU NTT、GPU naive、GPU Barrett 和 GPU Montgomery 四个版本的正确性验证。扩展测试覆盖 `poly_n = 1, 2, 3, 7, 16, 31, 127, 513, 1000, 4097`，测试数据包括全 0、全 1、递增序列、接近模数 `MOD` 的序列、固定随机种子和稀疏多项式。小规模测试同时对比 O(n^2) 朴素卷积、CPU NTT 和三个 GPU 版本，大规模测试对比 CPU NTT 和三个 GPU 版本。最终 `--verify-only` 输出 `ALL VERIFY TESTS PASSED`，说明当前实现中的 NTT 变换、逆变换、点乘以及三种模乘路径均能得到与 CPU reference 一致的结果。

从性能结果看，小规模输入下 GPU 版本并不占优。例如 `transform_n = 1024` 和 `4096` 时，GPU end-to-end 时间仍然高于 CPU baseline。这主要是因为小规模 NTT 的 butterfly 数量有限，GPU 并行度尚未充分展开，而一次 GPU 卷积仍然需要承担 kernel launch、cudaMalloc/cudaFree、H2D/D2H 拷贝和同步等固定开销。因此在小规模场景中，kernel-only 计算时间虽然较短，但端到端时间容易被固定开销主导。本实验环境中 CUDA runtime version 为 11070，即 CUDA 11.7；CUDA driver version 为 12060，即 CUDA 12.6。

随着 `transform_n` 增大，GPU 版本的并行优势逐渐显现。从 `2^14` 开始，GPU naive 和 Barrett 的 end-to-end 时间已经低于 CPU baseline；在 `2^20` 和 `2^22` 规模下，GPU 版本获得了较明显的加速。以 `transform_n = 2^22` 为例，CPU NTT mean 为 `561.601700 ms`，GPU naive total mean 为 `28.393265 ms`，GPU Montgomery total mean 为 `24.196315 ms`。这说明当 butterfly 数量足够大时，每个 CUDA thread 负责一个 butterfly 的映射方式能够有效利用 GPU 的大规模并行计算能力。

三个 GPU 版本之间的表现并非在所有规模上单调一致。本文统一使用 `naive_kernel / Barrett_kernel` 和 `naive_kernel / Montgomery_kernel` 表示优化版相对 naive 的 kernel 加速比；加速比大于 1 表示优化版 kernel 更快，加速比小于 1 表示优化版 kernel 更慢。naive 版本实现简单，部分规模下已经有较好表现；Barrett 版本在中小规模部分测试中略优于 naive，但在大规模上并未稳定领先；Montgomery 版本在 `2^18`、`2^20` 和 `2^22` 等较大规模下表现更好，但在小规模中会受到域转换和额外 kernel 开销影响。在当前实现中，由于 twiddle 计算、域转换或额外 kernel 开销，Barrett/Montgomery 并未在所有规模上超过 naive，但它们为后续减少整数取模开销提供了优化方向。后续若加入 twiddle 预处理、shared memory 小 stage 合并或 CUDA Graph，有望进一步降低重复计算和 kernel launch 开销。本机 nvcc 为 CUDA 11.7，无法直接使用 RTX 4060 对应的 sm_89 编译目标，因此 build.bat 自动回退到通用 CUDA 编译命令；本实验结果可以反映当前本地环境下的真实性能，但后续若升级 CUDA 11.8/12.x 并使用 sm_89 重新编译，可能进一步改善性能。
