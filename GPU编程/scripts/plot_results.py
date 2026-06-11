from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "bench.csv"
FIG_DIR = ROOT / "results" / "figures"


def setup_chinese_font():
    """设置常见中文字体，保证 Windows 本地绘图标题和图例尽量正常显示。"""
    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False


def pivot_metric(df, metric):
    """按 transform_n 展开不同算法的指标列，便于折线图绘制。"""
    table = df.pivot(index="transform_n", columns="algo", values=metric)
    return table.sort_index()


def plot_lines(table, title, ylabel, filename, log_y=False):
    """绘制多算法折线图。"""
    plt.figure(figsize=(9, 5.2))
    for algo in ["naive", "barrett", "montgomery"]:
        if algo in table.columns:
            plt.plot(table.index, table[algo], marker="o", linewidth=2, label=algo)
    plt.xlabel("transform_n")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend()
    plt.xscale("log", base=2)
    if log_y:
        plt.yscale("log")
    plt.tight_layout()
    plt.savefig(FIG_DIR / filename, dpi=180)
    plt.close()


def main():
    setup_chinese_font()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(CSV_PATH)
    df["transfer_ratio"] = (
        df["gpu_h2d_ms_mean"] + df["gpu_d2h_ms_mean"]
    ) / df["gpu_total_ms_mean"]

    plot_lines(
        pivot_metric(df, "gpu_total_ms_mean"),
        "GPU 端到端总时间对比（y 轴为 log scale）",
        "时间 / ms",
        "fig_total_time.png",
        log_y=True,
    )

    plot_lines(
        pivot_metric(df, "gpu_kernel_ms_mean"),
        "GPU kernel-only 时间对比（y 轴为 log scale）",
        "时间 / ms",
        "fig_kernel_time.png",
        log_y=True,
    )

    plot_lines(
        pivot_metric(df, "speedup_vs_cpu"),
        "GPU 相对 CPU 的端到端加速比",
        "加速比",
        "fig_speedup_vs_cpu.png",
        log_y=False,
    )

    plot_lines(
        pivot_metric(df, "transfer_ratio"),
        "H2D/D2H 数据传输占端到端时间比例",
        "传输占比",
        "fig_transfer_ratio.png",
        log_y=False,
    )

    plot_lines(
        pivot_metric(df, "throughput_mbutterfly_per_s"),
        "GPU butterfly 吞吐率对比",
        "吞吐率 / MButterfly/s",
        "fig_throughput.png",
        log_y=False,
    )

    print(f"figures saved to: {FIG_DIR}")


if __name__ == "__main__":
    main()
