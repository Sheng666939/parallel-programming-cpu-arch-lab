#include "ntt_cuda.cuh"

#include "common.hpp"
#include "mod_arith.cuh"
#include "timer.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <iostream>

namespace gpu_ntt {
namespace {

constexpr int BLOCK_SIZE = 256;

// 当前函数内使用的 CUDA 错误处理：记录错误后跳转到统一清理逻辑。
#define CUDA_CHECK_LOCAL(call)                                                       \
    do {                                                                            \
        cudaError_t err__ = (call);                                                  \
        if (err__ != cudaSuccess) {                                                  \
            std::cerr << "CUDA error: " << cudaGetErrorString(err__)                 \
                      << ", location " << __FILE__ << ":" << __LINE__ << "\n";      \
            ok = false;                                                             \
            goto cleanup;                                                           \
        }                                                                           \
    } while (0)

#define CUDA_KERNEL_CHECK_LOCAL() CUDA_CHECK_LOCAL(cudaGetLastError())

// 计算下标 i 在 log_n 位二进制中的 bit-reversal 结果。
__device__ __forceinline__ int bit_reverse_index(int i, int log_n) {
    int r = 0;
    for (int k = 0; k < log_n; ++k) {
        r = (r << 1) | (i & 1);
        i >>= 1;
    }
    return r;
}

// 线程映射：每个线程负责一个下标 i；当 i < rev(i) 时交换两个元素。
__global__ void bit_reverse_kernel(uint32_t* a, int n, int log_n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    int j = bit_reverse_index(i, log_n);
    if (i < j) {
        uint32_t tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
    }
}

// 线程映射：每个线程负责一个 butterfly。
// tid -> group = tid / mid, offset = tid % mid，
// pos1 = group * 2 * mid + offset, pos2 = pos1 + mid。
__global__ void ntt_stage_kernel(uint32_t* a,
                                 int n,
                                 int mid,
                                 uint32_t stage_root,
                                 int algo) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total_butterflies = n / 2;
    if (tid >= total_butterflies) {
        return;
    }

    int group = tid / mid;
    int offset = tid % mid;
    int pos1 = group * 2 * mid + offset;
    int pos2 = pos1 + mid;

    uint32_t w = 1u;
    if (algo == static_cast<int>(GpuNttAlgo::Montgomery)) {
        uint32_t root_mont = to_mont_value(stage_root);
        w = pow_mont_device(root_mont, static_cast<uint32_t>(offset));
    } else {
        w = pow_mod_device(stage_root, static_cast<uint32_t>(offset), algo);
    }

    uint32_t u = a[pos1];
    uint32_t v = 0;
    if (algo == static_cast<int>(GpuNttAlgo::Barrett)) {
        v = mul_mod_barrett(a[pos2], w);
    } else if (algo == static_cast<int>(GpuNttAlgo::Montgomery)) {
        v = mont_mul(a[pos2], w);
    } else {
        v = mul_mod_naive(a[pos2], w);
    }

    a[pos1] = add_mod(u, v);
    a[pos2] = sub_mod(u, v);
}

// 线程映射：每个线程负责一个系数的点值乘法。
__global__ void pointwise_mul_kernel(uint32_t* a,
                                     const uint32_t* b,
                                     int n,
                                     int algo) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    if (algo == static_cast<int>(GpuNttAlgo::Barrett)) {
        a[i] = mul_mod_barrett(a[i], b[i]);
    } else if (algo == static_cast<int>(GpuNttAlgo::Montgomery)) {
        a[i] = mont_mul(a[i], b[i]);
    } else {
        a[i] = mul_mod_naive(a[i], b[i]);
    }
}

// 线程映射：每个线程负责一个系数，将普通整数域转换到 Montgomery 域。
__global__ void to_mont_kernel(uint32_t* a, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        a[i] = to_mont_value(a[i]);
    }
}

// 线程映射：每个线程负责一个系数，将 Montgomery 域结果转换回普通整数域。
__global__ void from_mont_kernel(uint32_t* a, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        a[i] = from_mont_value(a[i]);
    }
}

// 线程映射：每个线程负责一个系数，在逆 NTT 后乘 n^{-1}。
__global__ void mul_scalar_kernel(uint32_t* a, int n, uint32_t scalar, int algo) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    if (algo == static_cast<int>(GpuNttAlgo::Barrett)) {
        a[i] = mul_mod_barrett(a[i], scalar);
    } else if (algo == static_cast<int>(GpuNttAlgo::Montgomery)) {
        uint32_t scalar_mont = to_mont_value(scalar);
        a[i] = mont_mul(a[i], scalar_mont);
    } else {
        a[i] = mul_mod_naive(a[i], scalar);
    }
}

int log2_exact(int n) {
    int log_n = 0;
    while ((1 << log_n) < n) {
        ++log_n;
    }
    return log_n;
}

void launch_ntt(uint32_t* d_a, int n, bool inverse, GpuNttAlgo algo, bool& ok) {
    const int log_n = log2_exact(n);
    const int blocks_n = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    const int algo_id = static_cast<int>(algo);
    bit_reverse_kernel<<<blocks_n, BLOCK_SIZE>>>(d_a, n, log_n);
    CUDA_KERNEL_CHECK_LOCAL();

    for (int len = 2; len <= n; len <<= 1) {
        const int mid = len >> 1;
        uint32_t stage_root = mod_pow(G, (MOD - 1) / len);
        if (inverse) {
            stage_root = mod_pow(stage_root, MOD - 2);
        }

        const int total_butterflies = n / 2;
        const int blocks = (total_butterflies + BLOCK_SIZE - 1) / BLOCK_SIZE;
        ntt_stage_kernel<<<blocks, BLOCK_SIZE>>>(d_a, n, mid, stage_root, algo_id);
        CUDA_KERNEL_CHECK_LOCAL();
    }

cleanup:
    return;
}

}  // namespace

const char* algo_name(GpuNttAlgo algo) {
    switch (algo) {
        case GpuNttAlgo::NaiveMod:
            return "naive";
        case GpuNttAlgo::Barrett:
            return "barrett";
        case GpuNttAlgo::Montgomery:
            return "montgomery";
        default:
            return "unknown";
    }
}

std::string get_cuda_device_name() {
    int device = 0;
    cudaError_t err = cudaGetDevice(&device);
    if (err != cudaSuccess) {
        return "Unknown CUDA Device";
    }
    cudaDeviceProp prop{};
    err = cudaGetDeviceProperties(&prop, device);
    if (err != cudaSuccess) {
        return "Unknown CUDA Device";
    }
    return prop.name;
}

int get_cuda_runtime_version() {
    int version = 0;
    cudaRuntimeGetVersion(&version);
    return version;
}

int get_cuda_driver_version() {
    int version = 0;
    cudaDriverGetVersion(&version);
    return version;
}

std::string get_cuda_compute_capability() {
    int device = 0;
    cudaError_t err = cudaGetDevice(&device);
    if (err != cudaSuccess) {
        return "unknown";
    }
    cudaDeviceProp prop{};
    err = cudaGetDeviceProperties(&prop, device);
    if (err != cudaSuccess) {
        return "unknown";
    }
    return std::to_string(prop.major) + "." + std::to_string(prop.minor);
}

int get_cuda_block_size() {
    return BLOCK_SIZE;
}

GpuNttResult convolution_cuda(const std::vector<uint32_t>& a,
                              const std::vector<uint32_t>& b,
                              std::vector<uint32_t>& out,
                              GpuNttAlgo algo) {
    GpuNttResult result{};
    HostTimer total_timer;
    bool ok = true;

    if (a.empty() || b.empty()) {
        out.clear();
        result.correct = true;
        return result;
    }

    const int result_size = static_cast<int>(a.size() + b.size() - 1);
    const int n = next_power_of_two(result_size);
    std::vector<uint32_t> ha(n, 0), hb(n, 0);
    std::copy(a.begin(), a.end(), ha.begin());
    std::copy(b.begin(), b.end(), hb.begin());
    out.assign(result_size, 0);

    uint32_t* d_a = nullptr;
    uint32_t* d_b = nullptr;
    cudaEvent_t h2d_start = nullptr;
    cudaEvent_t h2d_stop = nullptr;
    cudaEvent_t kernel_start = nullptr;
    cudaEvent_t kernel_stop = nullptr;
    cudaEvent_t d2h_start = nullptr;
    cudaEvent_t d2h_stop = nullptr;
    const int blocks_n = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    const int algo_id = static_cast<int>(algo);
    const uint32_t inv_n = mod_pow(static_cast<uint32_t>(n), MOD - 2);

    CUDA_CHECK_LOCAL(cudaEventCreate(&h2d_start));
    CUDA_CHECK_LOCAL(cudaEventCreate(&h2d_stop));
    CUDA_CHECK_LOCAL(cudaEventCreate(&kernel_start));
    CUDA_CHECK_LOCAL(cudaEventCreate(&kernel_stop));
    CUDA_CHECK_LOCAL(cudaEventCreate(&d2h_start));
    CUDA_CHECK_LOCAL(cudaEventCreate(&d2h_stop));

    CUDA_CHECK_LOCAL(cudaMalloc(&d_a, sizeof(uint32_t) * n));
    CUDA_CHECK_LOCAL(cudaMalloc(&d_b, sizeof(uint32_t) * n));

    CUDA_CHECK_LOCAL(cudaEventRecord(h2d_start));
    CUDA_CHECK_LOCAL(cudaMemcpy(d_a, ha.data(), sizeof(uint32_t) * n, cudaMemcpyHostToDevice));
    CUDA_CHECK_LOCAL(cudaMemcpy(d_b, hb.data(), sizeof(uint32_t) * n, cudaMemcpyHostToDevice));
    CUDA_CHECK_LOCAL(cudaEventRecord(h2d_stop));
    CUDA_CHECK_LOCAL(cudaEventSynchronize(h2d_stop));
    CUDA_CHECK_LOCAL(cudaEventElapsedTime(&result.h2d_ms, h2d_start, h2d_stop));

    CUDA_CHECK_LOCAL(cudaEventRecord(kernel_start));
    if (algo == GpuNttAlgo::Montgomery) {
        to_mont_kernel<<<blocks_n, BLOCK_SIZE>>>(d_a, n);
        CUDA_KERNEL_CHECK_LOCAL();
        to_mont_kernel<<<blocks_n, BLOCK_SIZE>>>(d_b, n);
        CUDA_KERNEL_CHECK_LOCAL();
    }

    launch_ntt(d_a, n, false, algo, ok);
    if (!ok) {
        goto cleanup;
    }
    launch_ntt(d_b, n, false, algo, ok);
    if (!ok) {
        goto cleanup;
    }

    pointwise_mul_kernel<<<blocks_n, BLOCK_SIZE>>>(d_a, d_b, n, algo_id);
    CUDA_KERNEL_CHECK_LOCAL();

    launch_ntt(d_a, n, true, algo, ok);
    if (!ok) {
        goto cleanup;
    }

    mul_scalar_kernel<<<blocks_n, BLOCK_SIZE>>>(d_a, n, inv_n, algo_id);
    CUDA_KERNEL_CHECK_LOCAL();

    if (algo == GpuNttAlgo::Montgomery) {
        from_mont_kernel<<<blocks_n, BLOCK_SIZE>>>(d_a, n);
        CUDA_KERNEL_CHECK_LOCAL();
    }
    CUDA_CHECK_LOCAL(cudaEventRecord(kernel_stop));
    CUDA_CHECK_LOCAL(cudaEventSynchronize(kernel_stop));
    CUDA_CHECK_LOCAL(cudaEventElapsedTime(&result.kernel_ms, kernel_start, kernel_stop));

    CUDA_CHECK_LOCAL(cudaEventRecord(d2h_start));
    CUDA_CHECK_LOCAL(cudaMemcpy(out.data(), d_a, sizeof(uint32_t) * result_size, cudaMemcpyDeviceToHost));
    CUDA_CHECK_LOCAL(cudaEventRecord(d2h_stop));
    CUDA_CHECK_LOCAL(cudaEventSynchronize(d2h_stop));
    CUDA_CHECK_LOCAL(cudaEventElapsedTime(&result.d2h_ms, d2h_start, d2h_stop));

cleanup:
    if (d_a != nullptr) {
        cudaFree(d_a);
    }
    if (d_b != nullptr) {
        cudaFree(d_b);
    }
    if (h2d_start != nullptr) {
        cudaEventDestroy(h2d_start);
    }
    if (h2d_stop != nullptr) {
        cudaEventDestroy(h2d_stop);
    }
    if (kernel_start != nullptr) {
        cudaEventDestroy(kernel_start);
    }
    if (kernel_stop != nullptr) {
        cudaEventDestroy(kernel_stop);
    }
    if (d2h_start != nullptr) {
        cudaEventDestroy(d2h_start);
    }
    if (d2h_stop != nullptr) {
        cudaEventDestroy(d2h_stop);
    }

    result.total_ms = total_timer.elapsed_ms();
    result.correct = ok;
    return result;
}

}  // namespace gpu_ntt
