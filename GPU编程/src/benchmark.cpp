#include "common.hpp"
#include "ntt_cpu.hpp"
#include "ntt_cuda.cuh"
#include "timer.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <numeric>
#include <string>
#include <vector>

namespace gpu_ntt {
namespace {

constexpr int CPU_NAIVE_CHECK_MAX_N = 4096;

struct BenchRow {
    std::string timestamp;
    std::string gpu_name;
    int cuda_runtime_version = 0;
    int cuda_driver_version = 0;
    std::string compute_capability;
    int transform_n = 0;
    int poly_n = 0;
    GpuNttAlgo algo = GpuNttAlgo::NaiveMod;
    int repeat = 0;
    int warmup = 0;
    int block_size = 0;
    double cpu_ms_mean = 0.0;
    double cpu_ms_median = 0.0;
    double gpu_h2d_ms = 0.0;
    double gpu_kernel_ms = 0.0;
    double gpu_d2h_ms = 0.0;
    double gpu_total_ms = 0.0;
    double gpu_total_ms_median = 0.0;
    double speedup_vs_cpu = 0.0;
    double speedup_vs_gpu_naive = 0.0;
    double throughput_mbutterfly_per_s = 0.0;
    bool correct = false;
};

std::vector<uint32_t> random_poly(int n, uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<uint32_t> dist(0, MOD - 1);
    std::vector<uint32_t> a(n);
    for (uint32_t& x : a) {
        x = dist(rng);
    }
    return a;
}

double ntt_butterflies(int transform_n) {
    const double log_n = std::log2(static_cast<double>(transform_n));
    return 3.0 * (static_cast<double>(transform_n) / 2.0) * log_n;
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const size_t mid = values.size() / 2;
    if (values.size() % 2 == 1) {
        return values[mid];
    }
    return (values[mid - 1] + values[mid]) / 2.0;
}

void write_csv(const std::string& path, const std::vector<BenchRow>& rows) {
    const bool exists = std::filesystem::exists(path);
    std::ofstream fout(path, std::ios::app);
    if (!exists) {
        fout << "timestamp,gpu_name,cuda_runtime_version,cuda_driver_version,compute_capability,"
             << "transform_n,poly_n,algo,repeat,warmup,block_size,"
             << "cpu_ms_mean,cpu_ms_median,gpu_h2d_ms_mean,gpu_kernel_ms_mean,"
             << "gpu_d2h_ms_mean,gpu_total_ms_mean,gpu_total_ms_median,"
             << "speedup_vs_cpu,speedup_vs_gpu_naive,throughput_mbutterfly_per_s,correct\n";
    }

    fout << std::fixed << std::setprecision(6);
    for (const BenchRow& row : rows) {
        fout << row.timestamp << ','
             << '"' << row.gpu_name << '"' << ','
             << row.cuda_runtime_version << ','
             << row.cuda_driver_version << ','
             << row.compute_capability << ','
             << row.transform_n << ','
             << row.poly_n << ','
             << algo_name(row.algo) << ','
             << row.repeat << ','
             << row.warmup << ','
             << row.block_size << ','
             << row.cpu_ms_mean << ','
             << row.cpu_ms_median << ','
             << row.gpu_h2d_ms << ','
             << row.gpu_kernel_ms << ','
             << row.gpu_d2h_ms << ','
             << row.gpu_total_ms << ','
             << row.gpu_total_ms_median << ','
             << row.speedup_vs_cpu << ','
             << row.speedup_vs_gpu_naive << ','
             << row.throughput_mbutterfly_per_s << ','
             << (row.correct ? "true" : "false") << '\n';
    }
}

BenchRow run_one_algo(const std::vector<uint32_t>& a,
                      const std::vector<uint32_t>& b,
                      const std::vector<uint32_t>& cpu_ref,
                      const BenchConfig& config,
                      int transform_n,
                      GpuNttAlgo algo,
                      const std::string& gpu_name,
                      int cuda_runtime_version,
                      int cuda_driver_version,
                      const std::string& compute_capability,
                      int block_size,
                      const std::string& timestamp) {
    std::vector<uint32_t> gpu_out;
    for (int i = 0; i < config.warmup; ++i) {
        (void)convolution_cuda(a, b, gpu_out, algo);
    }

    double h2d = 0.0;
    double kernel = 0.0;
    double d2h = 0.0;
    double total = 0.0;
    std::vector<double> total_samples;
    bool cuda_ok = true;
    for (int i = 0; i < config.repeat; ++i) {
        GpuNttResult r = convolution_cuda(a, b, gpu_out, algo);
        cuda_ok = cuda_ok && r.correct;
        h2d += r.h2d_ms;
        kernel += r.kernel_ms;
        d2h += r.d2h_ms;
        total += r.total_ms;
        total_samples.push_back(r.total_ms);
    }

    const double denom = std::max(1, config.repeat);
    BenchRow row;
    row.timestamp = timestamp;
    row.gpu_name = gpu_name;
    row.cuda_runtime_version = cuda_runtime_version;
    row.cuda_driver_version = cuda_driver_version;
    row.compute_capability = compute_capability;
    row.transform_n = transform_n;
    row.poly_n = static_cast<int>(a.size());
    row.algo = algo;
    row.repeat = config.repeat;
    row.warmup = config.warmup;
    row.block_size = block_size;
    row.gpu_h2d_ms = h2d / denom;
    row.gpu_kernel_ms = kernel / denom;
    row.gpu_d2h_ms = d2h / denom;
    row.gpu_total_ms = total / denom;
    row.gpu_total_ms_median = median(total_samples);
    row.correct = cuda_ok && compare_vectors(cpu_ref, gpu_out, algo_name(algo), true);
    if (row.gpu_kernel_ms > 0.0) {
        row.throughput_mbutterfly_per_s = ntt_butterflies(transform_n) / (row.gpu_kernel_ms * 1000.0);
    }
    return row;
}

struct VerifyCase {
    std::string name;
    std::vector<uint32_t> a;
    std::vector<uint32_t> b;
};

std::vector<uint32_t> make_sparse_poly(int n) {
    std::vector<uint32_t> a(n, 0);
    if (n > 0) {
        a[0] = 1;
        a[n / 2] = 7;
        a[n - 1] = MOD - 3;
    }
    for (int i = 5; i < n; i += 97) {
        a[i] = static_cast<uint32_t>((i * 17ull + 11ull) % MOD);
    }
    return a;
}

std::vector<VerifyCase> make_verify_cases(int n) {
    std::vector<uint32_t> zero(n, 0);
    std::vector<uint32_t> one(n, 1);
    std::vector<uint32_t> inc(n);
    std::vector<uint32_t> near_mod(n);
    for (int i = 0; i < n; ++i) {
        inc[i] = static_cast<uint32_t>(i % MOD);
        near_mod[i] = static_cast<uint32_t>((MOD - 1ull - static_cast<uint32_t>(i % MOD)) % MOD);
    }
    std::vector<uint32_t> rnd = random_poly(n, 202606u + static_cast<uint32_t>(n));
    std::vector<uint32_t> sparse = make_sparse_poly(n);

    return {
        {"all_zero", zero, zero},
        {"all_one", one, one},
        {"increasing", inc, inc},
        {"near_mod", near_mod, inc},
        {"random_seed_202606", rnd, random_poly(n, 202606u)},
        {"sparse", sparse, inc}
    };
}

bool report_first_mismatch(const std::vector<uint32_t>& expected,
                           const std::vector<uint32_t>& actual,
                           const std::string& case_name,
                           int poly_n,
                           const std::string& algo) {
    if (expected.size() != actual.size()) {
        std::cerr << "VERIFY FAIL case=" << case_name
                  << " poly_n=" << poly_n
                  << " algo=" << algo
                  << " size expected=" << expected.size()
                  << " actual=" << actual.size() << "\n";
        return false;
    }
    for (size_t i = 0; i < expected.size(); ++i) {
        if (expected[i] != actual[i]) {
            std::cerr << "VERIFY FAIL case=" << case_name
                      << " poly_n=" << poly_n
                      << " algo=" << algo
                      << " index=" << i
                      << " CPU=" << expected[i]
                      << " GPU=" << actual[i] << "\n";
            return false;
        }
    }
    return true;
}

}  // namespace

int run_verify_only() {
    const std::vector<int> sizes = {1, 2, 3, 7, 16, 31, 127, 513, 1000, 4097};
    bool all_ok = true;
    std::cout << "VERIFY GPU NTT convolution\n";
    std::cout << "GPU: " << get_cuda_device_name()
              << ", runtime=" << get_cuda_runtime_version()
              << ", driver=" << get_cuda_driver_version()
              << ", cc=" << get_cuda_compute_capability() << "\n";

    for (int poly_n : sizes) {
        for (const VerifyCase& test : make_verify_cases(poly_n)) {
            std::vector<uint32_t> cpu_ref = convolution_cpu(test.a, test.b);
            if (poly_n <= 1000) {
                std::vector<uint32_t> naive_ref = convolution_naive(test.a, test.b);
                all_ok = report_first_mismatch(naive_ref, cpu_ref, test.name, poly_n, "cpu_ntt") && all_ok;
            }

            for (GpuNttAlgo algo : {GpuNttAlgo::NaiveMod, GpuNttAlgo::Barrett, GpuNttAlgo::Montgomery}) {
                std::vector<uint32_t> gpu_out;
                GpuNttResult r = convolution_cuda(test.a, test.b, gpu_out, algo);
                const bool ok = r.correct &&
                    report_first_mismatch(cpu_ref, gpu_out, test.name, poly_n, algo_name(algo));
                all_ok = ok && all_ok;
            }
        }
        std::cout << "poly_n=" << poly_n << " verify done\n";
    }

    if (all_ok) {
        std::cout << "ALL VERIFY TESTS PASSED\n";
        return 0;
    }
    std::cout << "VERIFY TESTS FAILED\n";
    return 1;
}

int run_benchmark(const BenchConfig& config) {
    std::filesystem::create_directories("results");
    const std::string csv_path = "results/bench.csv";
    if (std::filesystem::exists(csv_path)) {
        std::filesystem::remove(csv_path);
    }
    const std::string gpu_name = get_cuda_device_name();
    const int cuda_runtime_version = get_cuda_runtime_version();
    const int cuda_driver_version = get_cuda_driver_version();
    const std::string compute_capability = get_cuda_compute_capability();
    const int block_size = get_cuda_block_size();
    const std::string timestamp = current_timestamp();

    std::cout << "GPU: " << gpu_name << "\n";
    std::cout << "CUDA runtime: " << cuda_runtime_version
              << ", driver: " << cuda_driver_version
              << ", compute capability: " << compute_capability
              << ", block size: " << block_size << "\n";
    std::cout << "CSV: " << csv_path << "\n";

    std::vector<BenchRow> all_rows;
    for (int log_n = config.min_log; log_n <= config.max_log; log_n += 2) {
        const int transform_n = 1 << log_n;
        const int poly_n = transform_n / 2;
        std::cout << "\n=== transform_n=2^" << log_n
                  << " (" << transform_n << "), poly_n=" << poly_n << " ===\n";

        std::vector<uint32_t> a = random_poly(poly_n, 1234u + log_n);
        std::vector<uint32_t> b = random_poly(poly_n, 5678u + log_n);

        HostTimer cpu_timer;
        std::vector<uint32_t> cpu_ref = convolution_cpu(a, b);
        const double cpu_ms = cpu_timer.elapsed_ms();
        std::cout << "CPU NTT: " << cpu_ms << " ms\n";

        if (transform_n <= CPU_NAIVE_CHECK_MAX_N) {
            std::vector<uint32_t> naive_ref = convolution_naive(a, b);
            const bool naive_ok = compare_vectors(naive_ref, cpu_ref, "CPU NTT vs O(n^2)", true);
            std::cout << "CPU vs O(n^2): " << (naive_ok ? "PASS" : "FAIL") << "\n";
        }

        std::vector<BenchRow> rows;
        rows.push_back(run_one_algo(a, b, cpu_ref, config, transform_n,
                                    GpuNttAlgo::NaiveMod, gpu_name,
                                    cuda_runtime_version, cuda_driver_version,
                                    compute_capability, block_size, timestamp));
        rows.push_back(run_one_algo(a, b, cpu_ref, config, transform_n,
                                    GpuNttAlgo::Barrett, gpu_name,
                                    cuda_runtime_version, cuda_driver_version,
                                    compute_capability, block_size, timestamp));
        rows.push_back(run_one_algo(a, b, cpu_ref, config, transform_n,
                                    GpuNttAlgo::Montgomery, gpu_name,
                                    cuda_runtime_version, cuda_driver_version,
                                    compute_capability, block_size, timestamp));

        const double naive_kernel = rows.front().gpu_kernel_ms;
        for (BenchRow& row : rows) {
            row.cpu_ms_mean = cpu_ms;
            row.cpu_ms_median = cpu_ms;
            if (row.gpu_total_ms > 0.0) {
                row.speedup_vs_cpu = cpu_ms / row.gpu_total_ms;
            }
            if (row.gpu_kernel_ms > 0.0 && naive_kernel > 0.0) {
                row.speedup_vs_gpu_naive = naive_kernel / row.gpu_kernel_ms;
            }
            std::cout << std::left << std::setw(11) << algo_name(row.algo)
                      << " correct=" << (row.correct ? "true " : "false")
                      << " h2d=" << row.gpu_h2d_ms
                      << " kernel=" << row.gpu_kernel_ms
                      << " d2h=" << row.gpu_d2h_ms
                      << " total=" << row.gpu_total_ms
                      << " ms\n";
        }

        all_rows.insert(all_rows.end(), rows.begin(), rows.end());
        write_csv(csv_path, rows);
    }

    std::cout << "\nbenchmark 完成，累计写入 " << all_rows.size() << " 行。\n";
    return 0;
}

}  // namespace gpu_ntt
