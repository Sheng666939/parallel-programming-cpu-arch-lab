cmake_minimum_required(VERSION 3.18)

project(gpu_ntt LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

add_executable(gpu_ntt
    src/main.cu
    src/ntt_cpu.cpp
    src/ntt_cuda.cu
    src/benchmark.cpp
)

target_include_directories(gpu_ntt PRIVATE include)

if (MSVC)
    target_compile_options(gpu_ntt PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/O2>)
endif()
