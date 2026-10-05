"""Fig 1 - Peak host memory of Ripes.exe against guest bytes written (Stage 1, measurement A).

Data: measure/results.md, RV32_ISS, Ripes v2.2.6-106-g5b8a616, median of 3 runs per point.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.environ.get(
    "KBFIG_HOME",
    r"C:\Users\itlab\AppData\Roaming\Claude\local-agent-mode-sessions\skills-plugin"
    r"\becb7c67-08ac-4617-82af-90aaf2fec4da\8b8fd7c4-f0fc-4d8d-ac53-cd809e148d29"
    r"\skills\kb-figure\scripts"))
from kbfig import Figure, load_table  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
df = load_table(os.path.join(HERE, "stage1_memory.csv"))

f = Figure(layout="single", palette="vivid", lang="en")
f.line(df, x="guest_MiB", y=["fill_peak_MiB", "control_peak_MiB"],
       labels=["fill (STRIDE = 4)", "control (STRIDE = 0)"],
       xlabel="Guest bytes written (MiB)",
       ylabel="Peak working set of Ripes.exe (MiB)",
       xlog=False, ylim=(0, 750), legend_loc="upper left")

# Least-squares fit through the fill points: slope is host MiB per guest MiB,
# which is the same number as host bytes per guest byte.
slope, intercept = np.polyfit(df["guest_MiB"], df["fill_peak_MiB"], 1)
xs = np.linspace(0, df["guest_MiB"].max(), 50)
f.ax.plot(xs, slope * xs + intercept, ls=":", lw=0.8, color="0.45", zorder=1)
f.ax.annotate(f"slope = {slope:.1f} host bytes / guest byte\n"
              f"intercept = {intercept:.0f} MiB (Ripes itself)",
              xy=(4, slope * 4 + intercept), xytext=(0.35, 510),
              fontsize=7, ha="left", va="center",
              arrowprops=dict(arrowstyle="-", lw=0.6, color="0.45"))

# Only the PNG is kept in the repo; it lands next to this script.
f.save("fig1_memory", outdir=HERE, formats=("png",), gnuplot=False, data=False)
print(f"fit: slope={slope:.2f}, intercept={intercept:.2f} MiB")
