"""A single plotting style for all notebooks: boxed axes, inward ticks."""

from __future__ import annotations

import matplotlib as mpl

__all__ = ["use_style", "panel_label"]


def use_style():
    """Apply the repository's Matplotlib style."""
    mpl.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 11,
            "legend.fontsize": 9.5,
            "legend.frameon": False,
            "axes.linewidth": 1.0,
            "axes.spines.top": True,
            "axes.spines.right": True,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.minor.visible": True,
            "ytick.minor.visible": True,
            "lines.linewidth": 1.6,
            "lines.markersize": 5.5,
            "errorbar.capsize": 2.5,
            "mathtext.fontset": "cm",
            "figure.constrained_layout.use": True,
        }
    )


def panel_label(ax, text, x=0.03, y=0.95):
    """Write a panel label such as ``(a)`` inside the top-left corner."""
    ax.text(x, y, text, transform=ax.transAxes, ha="left", va="top",
            fontweight="bold")
