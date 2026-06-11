#include "ntt_cpu.hpp"

#include "common.hpp"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace gpu_ntt {

uint32_t mod_pow(uint32_t base, uint64_t exp) {
    uint64_t result = 1;
    uint64_t cur = base;
    while (exp > 0) {
        if (exp & 1ull) {
            result = result * cur % MOD;
        }
        cur = cur * cur % MOD;
        exp >>= 1ull;
    }
    return static_cast<uint32_t>(result);
}

int next_power_of_two(int x) {
    int n = 1;
    while (n < x) {
        n <<= 1;
    }
    return n;
}

std::string current_timestamp() {
    const auto now = std::chrono::system_clock::now();
    const auto time = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &time);
#else
    localtime_r(&time, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

bool compare_vectors(const std::vector<uint32_t>& expected,
                     const std::vector<uint32_t>& actual,
                     const std::string& name,
                     bool verbose) {
    if (expected.size() != actual.size()) {
        if (verbose) {
            std::cerr << "[" << name << "] 长度不一致：expected=" << expected.size()
                      << ", actual=" << actual.size() << "\n";
        }
        return false;
    }
    for (size_t i = 0; i < expected.size(); ++i) {
        if (expected[i] != actual[i]) {
            if (verbose) {
                std::cerr << "[" << name << "] 结果错误：index=" << i
                          << ", CPU=" << expected[i]
                          << ", GPU=" << actual[i] << "\n";
            }
            return false;
        }
    }
    return true;
}

void ntt_cpu(std::vector<uint32_t>& a, bool invert) {
    const int n = static_cast<int>(a.size());

    // bit-reversal 重排：将自然顺序转换为迭代 Cooley-Tukey 所需的顺序。
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    // len 表示当前合并段长度，mid=len/2 是每个 butterfly 的跨度。
    for (int len = 2; len <= n; len <<= 1) {
        uint32_t wlen = mod_pow(G, (MOD - 1) / len);
        if (invert) {
            wlen = mod_pow(wlen, MOD - 2);
        }

        for (int i = 0; i < n; i += len) {
            uint32_t w = 1;
            const int mid = len >> 1;
            for (int j = 0; j < mid; ++j) {
                const uint32_t u = a[i + j];
                const uint32_t v = static_cast<uint32_t>(
                    static_cast<uint64_t>(a[i + j + mid]) * w % MOD);
                a[i + j] = u + v < MOD ? u + v : u + v - MOD;
                a[i + j + mid] = u >= v ? u - v : u + MOD - v;
                w = static_cast<uint32_t>(static_cast<uint64_t>(w) * wlen % MOD);
            }
        }
    }

    // 逆变换最后统一乘 n^{-1}，得到普通系数表示。
    if (invert) {
        const uint32_t inv_n = mod_pow(static_cast<uint32_t>(n), MOD - 2);
        for (uint32_t& x : a) {
            x = static_cast<uint32_t>(static_cast<uint64_t>(x) * inv_n % MOD);
        }
    }
}

std::vector<uint32_t> convolution_cpu(const std::vector<uint32_t>& a,
                                      const std::vector<uint32_t>& b) {
    if (a.empty() || b.empty()) {
        return {};
    }

    const int result_size = static_cast<int>(a.size() + b.size() - 1);
    const int n = next_power_of_two(result_size);
    std::vector<uint32_t> fa(n, 0), fb(n, 0);
    std::copy(a.begin(), a.end(), fa.begin());
    std::copy(b.begin(), b.end(), fb.begin());

    ntt_cpu(fa, false);
    ntt_cpu(fb, false);
    for (int i = 0; i < n; ++i) {
        fa[i] = static_cast<uint32_t>(static_cast<uint64_t>(fa[i]) * fb[i] % MOD);
    }
    ntt_cpu(fa, true);
    fa.resize(result_size);
    return fa;
}

std::vector<uint32_t> convolution_naive(const std::vector<uint32_t>& a,
                                        const std::vector<uint32_t>& b) {
    if (a.empty() || b.empty()) {
        return {};
    }
    std::vector<uint32_t> c(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            const uint64_t value = c[i + j] + static_cast<uint64_t>(a[i]) * b[j];
            c[i + j] = static_cast<uint32_t>(value % MOD);
        }
    }
    return c;
}

}  // namespace gpu_ntt
