#pragma once

#include <cstdint>
#include <vector>

namespace gpu_ntt {

// CPU 串行迭代 NTT。输入长度必须是 2 的幂。
void ntt_cpu(std::vector<uint32_t>& a, bool invert);

// 基于 CPU NTT 的多项式卷积 baseline。
std::vector<uint32_t> convolution_cpu(const std::vector<uint32_t>& a,
                                      const std::vector<uint32_t>& b);

// 小规模正确性验证用的 O(n^2) 朴素卷积。
std::vector<uint32_t> convolution_naive(const std::vector<uint32_t>& a,
                                        const std::vector<uint32_t>& b);

}  // namespace gpu_ntt
