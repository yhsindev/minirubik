"""Fig 2 - Simulation rate of RV32_ISS for the fill and control loops (Stage 1, measurement B).

Data: measure/results.md, Ripes v2.2.6-106-g5b8a616, IPS = instructions retired / --exectime,
median of 3 runs per point.
"""
import os
import sys

sys.path.insert(0, os.environ.get(
    "KBFIG_HOME",
    r"C:\Users\itlab\AppData\Roaming\Claude\local-agent-mode-sessions\skills-plugin"
    r"\becb7c67-08ac-4617-82af-90aaf2fec4da\8b8fd7c4-f0fc-4d8d-ac53-cd809e148d29"
    r"\skills\kb-figure\scripts"))
from kbfig import Figure, load_table  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
df = load_table(os.path.join(HERE, "stage1_ips.csv"))

f = Figure(layout="single", palette="vivid", lang="en")
f.line(df, x="guest_MiB", y=["control_MIPS", "fill_MIPS"],
       labels=["control (STRIDE = 0)", "fill (STRIDE = 4)"],
       xlabel="Guest bytes written (MiB)",
       ylabel="Simulation rate (million instr. / s)",
       ylim=(0, 14), legend_loc="center left")

ratio = (df["control_MIPS"] / df["fill_MIPS"])

# Only the PNG is kept in the repo; it lands next to this script.
f.save("fig2_ips", outdir=HERE, formats=("png",), gnuplot=False, data=False)
print("control/fill ratios (for the caption):", [round(r, 2) for r in ratio])
