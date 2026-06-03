# MPI NTT Polynomial Multiplication Lab

本目录用于保存并行程序设计 MPI 编程实验材料。实验选题为 NTT 多项式乘法优化，重点记录 MPI CRT 并行、普通输入 task-level MPI、CRT residue 收集方式对比，以及 MPI + rank 内线程混合并行版本。

## 实验信息

- 实验名称：MPI 编程实验
- 选题：NTT 多项式乘法 MPI 优化
- 平台：课程 OpenEuler MPI 环境
- 主要编译器：`mpic++`

## 主要版本

- `before_main.cc`：早期整理版本。
- `main_mpi_step1.cc`、`main_mpi_step1(1).cc`、`main_mpi_step2.cc`：MPI CRT 改造过程中的阶段版本。
- `main.cc`：稳定版本，支持普通输入 task-level MPI 与 input 4 的 CRT reduce-sum / point-to-point 收集切换。
- `main_v3_hybrid.cc`：Version 3，支持 MPI CRT + point-to-point collection + intra-rank threaded CRT NTT。

## 宏开关

- `THREAD_COUNT`：每个 MPI rank 使用的线程数。
- `ORDINARY_MPI_TASK=0`：input 0~3 使用 rank 0 Pthread whole-pipeline baseline。
- `ORDINARY_MPI_TASK=1`：input 0~3 使用 task-level MPI NTT，rank 0 forward A，rank 1 forward B。
- `CRT_COLLECT_METHOD=0`：input 4 使用 `MPI_Reduce + MPI_SUM` 收集 CRT residue。
- `CRT_COLLECT_METHOD=1`：input 4 使用 `MPI_Send / MPI_Recv` point-to-point 收集 CRT residue。
- `CRT_INTRA_THREAD=0`：Version 3 中 input 4 每个 CRT 小模数使用原单线程 NTT。
- `CRT_INTRA_THREAD=1`：Version 3 中 input 4 每个 CRT 小模数在 rank 内使用 OpenMP 多线程 NTT。

## 推荐编译命令

稳定 point-to-point 单线程版本：

```bash
mpic++ main_v3_hybrid.cc -O2 -fopenmp -pthread -DTHREAD_COUNT=1 -DORDINARY_MPI_TASK=0 -DCRT_COLLECT_METHOD=1 -DCRT_INTRA_THREAD=0 -o main_v3_ptp_t1
```

Hybrid 2 threads：

```bash
mpic++ main_v3_hybrid.cc -O2 -fopenmp -pthread -DTHREAD_COUNT=2 -DORDINARY_MPI_TASK=0 -DCRT_COLLECT_METHOD=1 -DCRT_INTRA_THREAD=1 -o main_v3_hybrid_t2
```

Hybrid 4 threads：

```bash
mpic++ main_v3_hybrid.cc -O2 -fopenmp -pthread -DTHREAD_COUNT=4 -DORDINARY_MPI_TASK=0 -DCRT_COLLECT_METHOD=1 -DCRT_INTRA_THREAD=1 -o main_v3_hybrid_t4
```

旧稳定版本仍可编译：

```bash
mpic++ main.cc -O2 -fopenmp -pthread -DTHREAD_COUNT=1 -DORDINARY_MPI_TASK=0 -DCRT_COLLECT_METHOD=1 -o main_ptp_t1
```

## 推荐运行组合

- `np=4, THREAD_COUNT=1`：point-to-point 单线程回退版。
- `np=2, THREAD_COUNT=2`：hybrid t2。
- `np=4, THREAD_COUNT=2`：当前最佳版本。
- `np=2, THREAD_COUNT=4`：hybrid t4。

资源约束按 `np * THREAD_COUNT <= nodes * ppn` 检查。例如：

- `np=4, THREAD_COUNT=1, nodes=1:ppn=4`，需要 4 核。
- `np=4, THREAD_COUNT=2, nodes=1:ppn=8`，需要 8 核。
- `np=2, THREAD_COUNT=4, nodes=1:ppn=8`，需要 8 核。
- `np=2, THREAD_COUNT=2, nodes=1:ppn=4`，需要 4 核。

## 当前最佳结果

input 4 当前最佳实验结果：

- `np=4`
- `THREAD_COUNT=2`
- `CRT_COLLECT_METHOD=1`
- `CRT_INTRA_THREAD=1`
- latency: `46.5961 ms`

结果整理见 `results/experiment_summary.md`。

## 注意事项

- 需要在课程 OpenEuler MPI 环境中通过 qsub 或课程指定方式运行。
- 不要直接在登录节点长时间运行性能测试。
- `/nttdata/` 输入和输出检查文件由课程环境提供，本仓库不提交课程数据文件。
- 可执行文件、qsub 临时输出和编译中间文件不提交。
