#pragma once

#include <cuda_runtime.h>
#include <cstdint>

namespace gpu_ntt {

static constexpr uint32_t MOD_DEVICE = 998244353u;
static constexpr uint32_t G_DEVICE = 3u;
static constexpr uint64_t BARRETT_MU = 18479187002ull;  // floor(2^64 / MOD)

static constexpr uint32_t MONT_N_PRIME = 998244351u;    // -MOD^{-1} mod 2^32
static constexpr uint32_t MONT_R = 301989884u;          // R mod MOD, R = 2^32
static constexpr uint32_t MONT_R2 = 932051910u;         // R^2 mod MOD

// 模加法：输入均在 [0, MOD) 范围内，返回 (a + b) mod MOD。
__device__ __forceinline__ uint32_t add_mod(uint32_t a, uint32_t b) {
    uint32_t s = a + b;
    if (s >= MOD_DEVICE || s < a) {
        s -= MOD_DEVICE;
    }
    return s;
}

// 模减法：输入均在 [0, MOD) 范围内，返回 (a - b) mod MOD。
__device__ __forceinline__ uint32_t sub_mod(uint32_t a, uint32_t b) {
    return (a >= b) ? (a - b) : (a + MOD_DEVICE - b);
}

// 朴素模乘：使用 64 位乘法后直接取模，作为 GPU baseline。
__device__ __forceinline__ uint32_t mul_mod_naive(uint32_t a, uint32_t b) {
    return static_cast<uint32_t>((static_cast<uint64_t>(a) * b) % MOD_DEVICE);
}

// Barrett 约简：用高 64 位乘法近似商，避免设备端整数除法。
__device__ __forceinline__ uint32_t barrett_reduce(uint64_t x) {
    uint64_t q = __umul64hi(x, BARRETT_MU);
    uint64_t r = x - q * MOD_DEVICE;
    if (r >= MOD_DEVICE) {
        r -= MOD_DEVICE;
    }
    if (r >= MOD_DEVICE) {
        r -= MOD_DEVICE;
    }
    return static_cast<uint32_t>(r);
}

// Barrett 模乘：适用于两个普通整数域元素的乘法。
__device__ __forceinline__ uint32_t mul_mod_barrett(uint32_t a, uint32_t b) {
    return barrett_reduce(static_cast<uint64_t>(a) * b);
}

// Montgomery 约简：输入 t < MOD * 2^32，返回 t * R^{-1} mod MOD。
__device__ __forceinline__ uint32_t mont_reduce(uint64_t t) {
    uint32_t m = static_cast<uint32_t>(t) * MONT_N_PRIME;
    uint64_t u = (t + static_cast<uint64_t>(m) * MOD_DEVICE) >> 32;
    if (u >= MOD_DEVICE) {
        u -= MOD_DEVICE;
    }
    return static_cast<uint32_t>(u);
}

// Montgomery 域乘法：若 a=xR, b=yR，则返回 xyR。
__device__ __forceinline__ uint32_t mont_mul(uint32_t a, uint32_t b) {
    return mont_reduce(static_cast<uint64_t>(a) * b);
}

// 将普通整数转换到 Montgomery 域：x -> xR。
__device__ __forceinline__ uint32_t to_mont_value(uint32_t x) {
    return mont_mul(x, MONT_R2);
}

// 将 Montgomery 域元素转换回普通整数域：xR -> x。
__device__ __forceinline__ uint32_t from_mont_value(uint32_t x) {
    return mont_reduce(x);
}

// 普通整数域快速幂，algo=0 使用朴素取模，algo=1 使用 Barrett。
__device__ __forceinline__ uint32_t pow_mod_device(uint32_t base, uint32_t exp, int algo) {
    uint32_t result = 1u;
    while (exp > 0u) {
        if (exp & 1u) {
            result = (algo == 1) ? mul_mod_barrett(result, base) : mul_mod_naive(result, base);
        }
        base = (algo == 1) ? mul_mod_barrett(base, base) : mul_mod_naive(base, base);
        exp >>= 1u;
    }
    return result;
}

// Montgomery 域快速幂，乘法单位元是 R mod MOD。
__device__ __forceinline__ uint32_t pow_mont_device(uint32_t base_mont, uint32_t exp) {
    uint32_t result = MONT_R;
    while (exp > 0u) {
        if (exp & 1u) {
            result = mont_mul(result, base_mont);
        }
        base_mont = mont_mul(base_mont, base_mont);
        exp >>= 1u;
    }
    return result;
}

}  // namespace gpu_ntt
