#pragma once

#include <chrono>

namespace gpu_ntt {

// 简单主机计时器，用于 CPU baseline 和 GPU end-to-end 总耗时统计。
class HostTimer {
public:
    HostTimer() { reset(); }

    void reset() {
        start_ = std::chrono::high_resolution_clock::now();
    }

    double elapsed_ms() const {
        const auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start_).count();
    }

private:
    std::chrono::high_resolution_clock::time_point start_;
};

}  // namespace gpu_ntt
