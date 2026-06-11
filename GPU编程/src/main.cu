#include "common.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

namespace gpu_ntt {
int run_benchmark(const BenchConfig& config);
int run_verify_only();
}

namespace {

void print_help() {
    std::cout << "Usage: gpu_ntt.exe [--verify-only] [--min-log 10] [--max-log 22] [--repeat 20] [--warmup 3]\n";
}

bool parse_int_arg(int argc, char** argv, int& i, int& value) {
    if (i + 1 >= argc) {
        return false;
    }
    value = std::atoi(argv[++i]);
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    gpu_ntt::BenchConfig config;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_help();
            return 0;
        }
        if (arg == "--verify-only") {
            config.verify_only = true;
        } else if (arg == "--min-log") {
            if (!parse_int_arg(argc, argv, i, config.min_log)) {
                std::cerr << "--min-log requires an integer argument\n";
                return 1;
            }
        } else if (arg == "--max-log") {
            if (!parse_int_arg(argc, argv, i, config.max_log)) {
                std::cerr << "--max-log requires an integer argument\n";
                return 1;
            }
        } else if (arg == "--repeat") {
            if (!parse_int_arg(argc, argv, i, config.repeat)) {
                std::cerr << "--repeat requires an integer argument\n";
                return 1;
            }
        } else if (arg == "--warmup") {
            if (!parse_int_arg(argc, argv, i, config.warmup)) {
                std::cerr << "--warmup requires an integer argument\n";
                return 1;
            }
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_help();
            return 1;
        }
    }

    if (config.min_log < 1 || config.max_log < config.min_log ||
        config.repeat < 1 || config.warmup < 0) {
        std::cerr << "Invalid arguments. Please check min/max log, repeat and warmup.\n";
        return 1;
    }

    if (config.verify_only) {
        return gpu_ntt::run_verify_only();
    }

    return gpu_ntt::run_benchmark(config);
}
