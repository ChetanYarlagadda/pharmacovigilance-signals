"""Shared matplotlib styling: colorblind-validated categorical palette, recessive chrome."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# Categorical slots, fixed order (validated for CVD separation on adjacent pairs).
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

# One colour per disproportionality method, used across every figure.
METHOD_COLOR = {"prr": SERIES[3], "ror": SERIES[1], "ic": SERIES[0], "ebgm": SERIES[2]}
METHOD_LABEL = {"prr": "PRR", "ror": "ROR", "ic": "BCPNN IC", "ebgm": "EBGM (GPS)"}
SCORE_LABEL = {"prr": "PRR", "prr025": "PRR lower 95%", "ror": "ROR", "ror025": "ROR lower 95%", "ic": "IC",
               "ic025": "IC lower 95%", "ebgm": "EBGM", "eb05": "EB05"}
SCORE_METHOD = {"prr": "prr", "prr025": "prr", "ror": "ror", "ror025": "ror", "ic": "ic", "ic025": "ic",
                "ebgm": "ebgm", "eb05": "ebgm"}
POSITIVE, NEGATIVE = SERIES[7], SERIES[0]  # positive / negative controls


def apply_style() -> None:
    """Set the shared rcParams for every report figure."""
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "figure.dpi": 150,
            "savefig.dpi": 150,
            "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"],
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.titlecolor": INK,
            "axes.labelcolor": INK_2,
            "axes.edgecolor": AXIS,
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.axisbelow": True,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelcolor": INK_2,
            "ytick.labelcolor": INK_2,
            "legend.frameon": False,
            "legend.labelcolor": INK_2,
            "lines.linewidth": 2,
            "lines.markersize": 6,
        }
    )


def save(fig: plt.Figure, path: Path) -> None:
    """Write a figure as a 150 dpi PNG and close it."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
