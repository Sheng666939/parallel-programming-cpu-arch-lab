#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace gpu_ntt {

// NTT 使用的常用 NTT 友好素数：MOD = 119 * 2^23 + 1。
inline constexpr uint32_t MOD = 998244353u;
inline constexpr uint32_t G = 3u;

// benchmark 默认参数，默认测试 2^10, 2^12, ..., 2^22。
struct BenchConfig {
    int min_log = 10;
    int max_log = 22;
    int repeat = 20;
    int warmup = 3;
    bool verify_only = false;
};

// 主机端快速幂，用于 CPU NTT 以及 GPU stage root 的预计算。
uint32_t mod_pow(uint32_t base, uint64_t exp);

// 计算不小于 x 的最小 2 的幂。
int next_power_of_two(int x);

// 返回当前时间字符串，写入 CSV 时用于标记实验批次。
std::string current_timestamp();

// 打印多项式结果不一致的位置，便于定位 GPU 算法错误。
bool compare_vectors(const std::vector<uint32_t>& expected,
                     const std::vector<uint32_t>& actual,
                     const std::string& name,
                     bool verbose);

}  // namespace gpu_ntt
