# MPI NTT 多项式乘法实验

这里整理的是我在 MPI 编程实验中完成的 NTT 多项式乘法优化代码。实验主要围绕两类输入展开：input 0~3 是普通小模数多项式乘法，input 4 是大模数，需要用多个 NTT 友好的小模数做 CRT，再合并回目标模数。

代码保留了原来的读取和检查接口，最终结果仍由 `fCheck` / `fCheckLarge` 判断。为了避免多进程重复输出，程序中的正确性结果和主要计时信息只由 rank 0 打印。

## 文件说明

- `before_main.cc`：MPI 改造前后保留下来的早期版本。
- `main_mpi_step1.cc`、`main_mpi_step1(1).cc`、`main_mpi_step2.cc`：MPI CRT 并行逐步修改过程中的阶段版本。
- `main.cc`：稳定版本，支持 input 0~3 的 task-level MPI 开关，以及 input 4 的 reduce-sum / point-to-point 收集方式切换。
- `main_v3_hybrid.cc`：完整保留版。在 input 4 上保留 MPI CRT 的 rank 间划分，同时给每个 rank 内部的单个 CRT 小模数 NTT 增加 OpenMP 线程并行。
- `main_v3_hybrid_final.cc`：精简后的最终版。它保留最终实验会用到的所有宏开关和检查流程，删去了不再被 `main()` 调用的旧本地 CRT 包装函数，便于阅读和提交。
- `scripts/qsub_mpi_ntt_template.sh`：qsub 运行模板，实际提交作业时可以按线程数和进程数调整。
- `results/experiment_summary.md`：整理后的关键实验结果。

## 实现路线

### 1. baseline 和 MPI CRT

input 0~3 先保留 rank 0 上的 Pthread whole-pipeline NTT baseline，这样普通小模数输入有一个稳定的正确性参考。

input 4 的模数比较大，不能直接用单个 NTT 友好模数完成，所以采用 4 个小模数分别做 NTT，再用 Garner/CRT 思路合并。MPI 层面的划分方式是：

```cpp
for (int modulusIndex = rank; modulusIndex < CRT_MOD_COUNT; modulusIndex += size)
```

也就是不同 rank 负责不同的 CRT 小模数。`np=4` 时每个 rank 基本负责一个小模数，`np=2` 时每个 rank 负责两个小模数，`np=1` 时退化为单进程处理全部小模数。

### 2. input 0~3 的 task-level MPI

普通小模数输入增加了一个任务级 MPI 版本：

- rank 0 做多项式 A 的 forward NTT；
- rank 1 做多项式 B 的 forward NTT；
- rank 1 把 B 的频域结果发给 rank 0；
- rank 0 完成点乘、inverse NTT 和检查。

这个版本主要用于展示 MPI 任务划分思路。由于后半段仍集中在 rank 0 上，性能收益不一定稳定，所以默认仍可以回到 baseline。

### 3. reduce-sum 和 point-to-point 收集

input 4 的 CRT residue 收集保留了两种方式：

- `CRT_COLLECT_METHOD=0`：使用 `MPI_Reduce + MPI_SUM`。每个 rank 把自己负责的小模数结果放到对应位置，其余位置为 0，最后在 rank 0 求和得到完整 residue 表。
- `CRT_COLLECT_METHOD=1`：使用 `MPI_Send / MPI_Recv`。每个 rank 只发送自己实际负责的小模数结果，rank 0 按 `modulusIndex` 接收并放到 `gatheredResidues[modulusIndex * resultLength + index]`。

point-to-point 方式避免了发送大量全 0 residue，input 4 上表现更好。

### 4. Version 3 hybrid

`main_v3_hybrid.cc` 是当前主要优化版本：

- MPI 层面仍按 CRT 小模数分配任务；
- 每个 rank 仍顺序处理自己负责的 `modulusIndex`，不会同时开多个 CRT 小模数的线程池；
- 单个 CRT 小模数内部用 OpenMP 并行清零/拷贝、NTT stage 内 butterfly、点乘和 residue 写回；
- bit reversal 暂时保留串行，避免 swap 写冲突；
- stage 之间依靠 OpenMP barrier 保证顺序。

这样资源使用比较清楚：总核数大致按 `np * THREAD_COUNT` 估算。

## 宏开关

- `THREAD_COUNT`：每个 MPI rank 内使用的线程数。
- `ORDINARY_MPI_TASK=0`：input 0~3 使用 rank 0 Pthread baseline。
- `ORDINARY_MPI_TASK=1`：input 0~3 使用 task-level MPI。
- `CRT_COLLECT_METHOD=0`：input 4 使用 reduce-sum 收集。
- `CRT_COLLECT_METHOD=1`：input 4 使用 point-to-point 收集。
- `CRT_INTRA_THREAD=0`：input 4 每个 CRT 小模数使用单线程 NTT。
- `CRT_INTRA_THREAD=1`：input 4 每个 CRT 小模数使用 rank 内 OpenMP 多线程 NTT。

## 编译命令

point-to-point 单线程版本：

```bash
mpic++ main_v3_hybrid_final.cc -O2 -fopenmp -pthread -DTHREAD_COUNT=1 -DORDINARY_MPI_TASK=0 -DCRT_COLLECT_METHOD=1 -DCRT_INTRA_THREAD=0 -o main_v3_ptp_t1
```

hybrid 2 线程版本：

```bash
mpic++ main_v3_hybrid_final.cc -O2 -fopenmp -pthread -DTHREAD_COUNT=2 -DORDINARY_MPI_TASK=0 -DCRT_COLLECT_METHOD=1 -DCRT_INTRA_THREAD=1 -o main_v3_hybrid_t2
```

hybrid 4 线程版本：

```bash
mpic++ main_v3_hybrid_final.cc -O2 -fopenmp -pthread -DTHREAD_COUNT=4 -DORDINARY_MPI_TASK=0 -DCRT_COLLECT_METHOD=1 -DCRT_INTRA_THREAD=1 -o main_v3_hybrid_t4
```

如果需要回到旧的稳定 MPI CRT 版本，也可以直接编译 `main.cc`：

```bash
mpic++ main.cc -O2 -fopenmp -pthread -DTHREAD_COUNT=1 -DORDINARY_MPI_TASK=0 -DCRT_COLLECT_METHOD=1 -o main_ptp_t1
```

## 推荐运行组合

- `np=4, THREAD_COUNT=1`：point-to-point 单线程回退版。
- `np=2, THREAD_COUNT=2`：hybrid t2，用 4 核。
- `np=4, THREAD_COUNT=2`：hybrid t2，用 8 核，目前记录中 input 4 最快。
- `np=2, THREAD_COUNT=4`：hybrid t4，用 8 核。

资源申请时按下面这个关系检查：

```text
np * THREAD_COUNT <= nodes * ppn
```

例如 `np=4, THREAD_COUNT=2` 时需要 8 个核心，可以使用 `nodes=1:ppn=8`。

## 当前最好结果

目前记录到的 input 4 最好结果为：

```text
np = 4
THREAD_COUNT = 2
CRT_COLLECT_METHOD = 1
CRT_INTRA_THREAD = 1
input 4 latency = 46.5961 ms
```

对应的是 `main_v3_hybrid_final.cc` 中的 MPI CRT + point-to-point + rank 内 2 线程版本。

## 运行注意

- 测试需要课程服务器上的 `/nttdata/` 数据和检查文件，本仓库没有提交这些数据。
- 性能测试应通过 qsub 作业运行，不建议在登录节点直接长时间跑。
- 可执行文件、qsub 临时输出和本地日志不放进仓库。
