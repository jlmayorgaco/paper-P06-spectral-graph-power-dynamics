"""Plot the P1 terminal-transfer conditioning diagnostics."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


CAMPAIGN = Path(__file__).resolve().parents[2]
raw = CAMPAIGN / "raw" / "gfl11" / "p1_terminal_transfer.csv"
out_dir = CAMPAIGN / "figures"
out_dir.mkdir(parents=True, exist_ok=True)
df = pd.read_csv(raw)

fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True, constrained_layout=True)
axes[0].loglog(df["frequency_hz"], df["condition"], color="#315b8a", lw=1.8)
axes[0].set_ylabel("cond(jωI − A)")
axes[0].grid(True, which="both", alpha=0.25)
axes[0].set_title("P1 frozen GFL11 terminal-transfer diagnostic")
axes[1].loglog(df["frequency_hz"], df["sigma_min"], color="#a33f2b", lw=1.8)
axes[1].set_xlabel("frequency [Hz]")
axes[1].set_ylabel("σ_min(jωI − A)")
axes[1].grid(True, which="both", alpha=0.25)
fig.savefig(out_dir / "P1_GFL_TRANSFER_CONDITIONING.png", dpi=180)
print(out_dir / "P1_GFL_TRANSFER_CONDITIONING.png")
