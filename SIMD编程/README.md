# Parallel Programming CPU Architecture Lab

本仓库用于保存并行程序设计课程 SIMD 编程实验相关材料。

## 实验信息

- 实验名称：SIMD 编程
- 选题：NTT 多项式乘法优化
- 平台：OpenEuler / ARM AArch64 / Kunpeng-920 / NEON
- 学生：罗晟皓
- 学号：2411566
- 专业：计算机科学与技术

## 文件说明

由于本次实验的主要代码编写、编译、运行和测试均在课程服务器上完成，实验过程中我对每一个重要版本的代码、汇编文件和指令统计结果都进行了归档保存。

本仓库中的若干 `.docs` 文件是从服务器实验记录中整理出的 Word 版本，用于记录不同阶段的代码、汇编片段、运行结果和分析过程。

其中主要包括：

- 基础标量 NTT 版本；
- NEON 蝴蝶加减优化版本；
- stageRoots 旋转因子预计算版本；
- NEON + stageRoots 稳定版本；
- Montgomery root-mul 标量优化版本；
- NEON Montgomery root-mul 实验版本；
- 第四个大模数 input 4 的 `long long + __int128` 探索版本；
- 各版本对应的汇编代码和指令统计结果。

最终稳定提交版本对应 `main.cc`，大模数探索版本对应 `main_bigmod_experiment.cc`。
