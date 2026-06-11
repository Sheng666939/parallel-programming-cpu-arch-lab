@echo off
setlocal

if not exist build mkdir build
if not exist bin mkdir bin

where cl >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
  if exist "%ProgramFiles%\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat" (
    call "%ProgramFiles%\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat" -arch=x64
  ) else if exist "%ProgramFiles%\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat" (
    call "%ProgramFiles%\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat" -arch=x64
  ) else (
    echo Cannot find cl.exe. Please run this script in x64 Native Tools Command Prompt for VS.
    exit /b 1
  )
)

echo [1/2] Try build with -arch=sm_89 ...
nvcc -O3 -std=c++17 -allow-unsupported-compiler -D_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH -Xcompiler /utf-8 -arch=sm_89 ^
  src/main.cu src/ntt_cpu.cpp src/ntt_cuda.cu src/benchmark.cpp ^
  -Iinclude ^
  -o bin/gpu_ntt.exe

if %ERRORLEVEL% EQU 0 (
  echo Build success: bin\gpu_ntt.exe
  exit /b 0
)

echo.
echo sm_89 build failed. Try generic CUDA build for older CUDA Toolkit ...
nvcc -O3 -std=c++17 -allow-unsupported-compiler -D_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH -Xcompiler /utf-8 ^
  src/main.cu src/ntt_cpu.cpp src/ntt_cuda.cu src/benchmark.cpp ^
  -Iinclude ^
  -o bin/gpu_ntt.exe

if %ERRORLEVEL% EQU 0 (
  echo Build success: bin\gpu_ntt.exe
  exit /b 0
)

echo Build failed. Please check CUDA Toolkit and MSVC environment.
exit /b 1
