# MPI NTT Experiment Summary

本文件记录 MPI NTT 实验中整理出的关键版本与结果，避免提交大量重复 qsub warning、临时 `.e/.o/.out/.log` 文件。

## Version Notes

| Version | Macro Settings | Description |
| --- | --- | --- |
| rank 0 baseline | `ORDINARY_MPI_TASK=0` | input 0~3 使用 rank 0 Pthread whole-pipeline baseline。 |
| task-level MPI | `ORDINARY_MPI_TASK=1` | input 0~3 中 rank 0 forward A，rank 1 forward B，rank 0 完成 pointwise 和 inverse。 |
| MPI CRT reduce-sum | `CRT_COLLECT_METHOD=0` | input 4 的 4 个 CRT 小模数按 rank 分配，使用 `MPI_Reduce + MPI_SUM` 收集 residue。 |
| MPI CRT point-to-point | `CRT_COLLECT_METHOD=1` | input 4 中每个 rank 只发送自己负责的小模数 residue，rank 0 使用 `MPI_Send / MPI_Recv` 收集。 |
| Version 3 hybrid | `CRT_COLLECT_METHOD=1`, `CRT_INTRA_THREAD=1` | input 4 使用 MPI CRT + point-to-point + rank 内 OpenMP threaded CRT NTT。最终提交版本为 `main_v3_hybrid_final.cc`。 |

## Best Known Result

| input | np | THREAD_COUNT | CRT_COLLECT_METHOD | CRT_INTRA_THREAD | latency (ms) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4 | 4 | 2 | 1 | 1 | 46.5961 |

## Recommended Runs

| Run | Resource Request | Reason |
| --- | --- | --- |
| `np=4, THREAD_COUNT=1` | `nodes=1:ppn=4` | `4 * 1 = 4 <= 4` |
| `np=2, THREAD_COUNT=2` | `nodes=1:ppn=4` | `2 * 2 = 4 <= 4` |
| `np=4, THREAD_COUNT=2` | `nodes=1:ppn=8` | `4 * 2 = 8 <= 8` |
| `np=2, THREAD_COUNT=4` | `nodes=1:ppn=8` | `2 * 4 = 8 <= 8` |

## Notes

- 所有最终正确性检查仍由 rank 0 调用 `fCheck` 或 `fCheckLarge`。
- `CRT_INTRA_THREAD=0` 可用于回退到旧单线程 CRT NTT 行为。
- qsub 临时输出和编译产物不进入仓库。
