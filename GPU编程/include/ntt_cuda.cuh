#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace gpu_ntt {

// GPU NTT 算法版本。Montgomery 版本作为高阶优化版本保留并参与测试。
enum class GpuNttAlgo {
    NaiveMod = 0,
    Barrett = 1,
    Montgomery = 2
};

struct GpuNttResult {
    bool correct = false;
    float h2d_ms = 0.0f;
    float kernel_ms = 0.0f;
    float d2h_ms = 0.0f;
    double total_ms = 0.0;
};

// 执行一次 GPU 多项式卷积，并分别返回 H2D、kernel-only、D2H 和 end-to-end 时间。
GpuNttResult convolution_cuda(const std::vector<uint32_t>& a,
                              const std::vector<uint32_t>& b,
                              std::vector<uint32_t>& out,
                              GpuNttAlgo algo);

// 查询当前 CUDA 设备名称，写入 benchmark CSV。
std::string get_cuda_device_name();

// 查询 CUDA runtime / driver 版本和设备 compute capability。
int get_cuda_runtime_version();
int get_cuda_driver_version();
std::string get_cuda_compute_capability();
int get_cuda_block_size();

const char* algo_name(GpuNttAlgo algo);

}  // namespace gpu_ntt
