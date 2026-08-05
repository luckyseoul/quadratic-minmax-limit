#!/usr/bin/env python3
"""
Render a 2-page summary of lim α_n = 1/2 for X / MO 413935.

Style: dark cards (#0f1419 / #1a2332 / accent #1d9bf0).
Output: x-cards/page1.jpg, page2.jpg, twopager.pdf (+ evidence/share mirrors)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "x-cards"
OUT.mkdir(exist_ok=True)

W_IN, H_IN = 8.5, 11.0
DPI = 220

BG = "#0f1419"
CARD = "#1a2332"
ACCENT = "#1d9bf0"
GREEN = "#00ba7c"
AMBER = "#ffad1f"
TEXT = "#e7e9ea"
MUTED = "#8b98a5"
WHITE = "#ffffff"
DIM = "#38444d"
GREEN_CARD = "#0a2a1f"
AMBER_CARD = "#2a2210"


def _fig():
    fig = plt.figure(figsize=(W_IN, H_IN), facecolor=BG, dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_facecolor(BG)
    return fig, ax


def _card(ax, x, y, w, h, fc=CARD):
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            linewidth=0,
            facecolor=fc,
            transform=ax.transAxes,
            clip_on=False,
        )
    )


def _bar(ax, y, color=ACCENT, x=0.06, w=0.88, h=0.006):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=color, edgecolor="none"))


def _t(ax, x, y, s, **kw):
    defaults = dict(color=TEXT, fontfamily="sans-serif", transform=ax.transAxes)
    defaults.update(kw)
    ax.text(x, y, s, **defaults)


def page1():
    """Prize-thread page 1: the answer, not a progress log."""
    fig, ax = _fig()

    # Header
    _t(ax, 0.06, 0.955, "MIN–MAX  ±1  QUADRATIC  FORM", fontsize=13, color=ACCENT, fontweight="bold")
    _t(
        ax,
        0.06,
        0.920,
        "MathOverflow 413935  ·  existence of the limit",
        fontsize=10,
        color=MUTED,
    )
    _bar(ax, 0.902)

    # Main theorem
    _card(ax, 0.05, 0.700, 0.90, 0.175, fc=GREEN_CARD)
    _t(ax, 0.08, 0.840, "THEOREM", fontsize=12, color=GREEN, fontweight="bold")
    _t(
        ax,
        0.50,
        0.775,
        r"$L \;=\; \lim_{n\to\infty}\alpha_n \;=\; \frac{1}{2}$",
        fontsize=22,
        color=WHITE,
        ha="center",
        va="center",
    )
    _t(
        ax,
        0.08,
        0.720,
        r"The limit exists and equals $\frac{1}{2}$.",
        fontsize=11,
        color=TEXT,
    )

    # The question
    _t(ax, 0.06, 0.655, "THE QUESTION", fontsize=11, color=AMBER, fontweight="bold")
    _t(
        ax,
        0.50,
        0.595,
        r"$L \;=\; \lim_{n\to\infty}\, n^{-3/2}\;"
        r"\min_{a_{ij}=\pm1}\;\max_{x_i=\pm1}"
        r"\left|\sum_{i<j} a_{ij}\,x_i x_j\right|$",
        fontsize=12.5,
        color=WHITE,
        ha="center",
        va="center",
    )
    _t(
        ax,
        0.06,
        0.545,
        r"Write $\alpha_n=m_n/n^{3/2}$.  Does $\lim\alpha_n$ exist?",
        fontsize=10.5,
        color=MUTED,
    )

    # Why 1/2
    _card(ax, 0.05, 0.250, 0.90, 0.260)
    _t(ax, 0.08, 0.475, "ARGUMENT", fontsize=11, color=GREEN, fontweight="bold")
    why = [
        r"1.  Sandwich:  $1/\pi \leq \liminf\alpha_n \leq \limsup\alpha_n \leq 1/2$.",
        r"2.  On the dense Paley family $n=p^2+1$,  $\Phi(C)=\frac{1}{2}n\sqrt{n-1}$.",
        r"     If $m_n\geq\Phi(C)-O(1)$, denseness forces $L=\frac{1}{2}$.",
        r"3.  Every gap-2 undercutter of $C$ is empty or non-descending;",
        r"     hence $m_n\geq\Phi(C)-2$ and $L=\frac{1}{2}$.",
    ]
    y = 0.430
    for line in why:
        _t(ax, 0.08, y, line, fontsize=10.5, color=TEXT)
        y -= 0.032

    # Scope
    _card(ax, 0.05, 0.095, 0.90, 0.120, fc=AMBER_CARD)
    _t(ax, 0.08, 0.180, "REMARK", fontsize=11, color=AMBER, fontweight="bold")
    _t(
        ax,
        0.08,
        0.130,
        "An independent residual / 16N spectral package remains open;\n"
        "it is not required for the limit.",
        fontsize=10.5,
        color=TEXT,
        linespacing=1.4,
    )

    _t(
        ax,
        0.06,
        0.050,
        "MathOverflow 413935  ·  github.com/luckyseoul/quadratic-minmax-limit",
        fontsize=9,
        color=MUTED,
    )
    _t(ax, 0.94, 0.050, "1 / 2", fontsize=10, color=MUTED, ha="right")
    _t(ax, 0.06, 0.022, "2026-08-05", fontsize=9, color=MUTED)

    path = OUT / "page1.jpg"
    fig.savefig(path, dpi=DPI, facecolor=BG, format="jpeg", pil_kwargs={"quality": 92})
    plt.close(fig)
    print("WROTE", path)


def page2():
    """Prize-thread page 2: what closed, no repo dump."""
    fig, ax = _fig()

    _t(ax, 0.06, 0.955, "OUTLINE OF THE PROOF", fontsize=14, color=ACCENT, fontweight="bold")
    _t(
        ax,
        0.06,
        0.920,
        r"bi-tight emptiness  ·  Type I  ·  deep freeness-fail  ·  denseness",
        fontsize=10,
        color=MUTED,
    )
    _bar(ax, 0.902)

    # Three pillars
    _card(ax, 0.05, 0.680, 0.90, 0.195)
    _t(ax, 0.08, 0.840, "①  BI-TIGHT COVERS", fontsize=11, color=GREEN, fontweight="bold")
    _t(
        ax,
        0.08,
        0.790,
        r"No size-$2p$ bi-tight Max$\pm$ cover of the Paley conference exists for primes $p\geq 5$.",
        fontsize=10.5,
        color=TEXT,
    )
    _t(
        ax,
        0.08,
        0.745,
        r"Multiplicity of $\lambda_{\max}$ and a lower bound on $\lambda_{\min}$ force the cycle",
        fontsize=10.5,
        color=TEXT,
    )
    _t(
        ax,
        0.08,
        0.705,
        r"eigenvalue below the bi-tight threshold (majorization; residual / 16N unused).",
        fontsize=10.5,
        color=TEXT,
    )

    _card(ax, 0.05, 0.470, 0.90, 0.185)
    _t(ax, 0.08, 0.620, "②  TYPE-I FREENESS FAILURE", fontsize=11, color=GREEN, fontweight="bold")
    _t(
        ax,
        0.08,
        0.570,
        r"The remaining Type-I gap-2 candidate ($k=3p-2$) dualizes to a linear system",
        fontsize=10.5,
        color=TEXT,
    )
    _t(
        ax,
        0.08,
        0.525,
        r"on edge variables.  Box and sum constraints are contradictory for every $p\geq 5$.",
        fontsize=10.5,
        color=TEXT,
    )
    _t(
        ax,
        0.08,
        0.485,
        "Hence that class admits no further descent.",
        fontsize=10.5,
        color=TEXT,
    )

    _card(ax, 0.05, 0.260, 0.90, 0.185)
    _t(ax, 0.08, 0.410, "③  DEEP FREENESS FAILURE", fontsize=11, color=GREEN, fontweight="bold")
    _t(
        ax,
        0.08,
        0.360,
        r"Deep covers ($k\geq 3p$) are either free, fail-equality empty, or dualize to a",
        fontsize=10.5,
        color=TEXT,
    )
    _t(
        ax,
        0.08,
        0.315,
        r"two-level system obstructed by the same box lower bound.",
        fontsize=10.5,
        color=TEXT,
    )
    _t(
        ax,
        0.08,
        0.275,
        r"With freeness and tight emptiness:  $m_n \geq \Phi(C)-2$.",
        fontsize=10.5,
        color=TEXT,
    )

    # Main claim
    _card(ax, 0.05, 0.125, 0.90, 0.105, fc=GREEN_CARD)
    _t(ax, 0.08, 0.195, "CONCLUSION", fontsize=11, color=GREEN, fontweight="bold")
    _t(
        ax,
        0.08,
        0.145,
        r"$L=\frac{1}{2}$.   Writeup and verification: github.com/luckyseoul/quadratic-minmax-limit",
        fontsize=11,
        color=WHITE,
    )

    _bar(ax, 0.095, color=DIM)
    _t(
        ax,
        0.06,
        0.050,
        "MathOverflow 413935  ·  github.com/luckyseoul/quadratic-minmax-limit",
        fontsize=9.5,
        color=MUTED,
    )
    _t(ax, 0.94, 0.050, "2 / 2", fontsize=10, color=MUTED, ha="right")
    _t(ax, 0.06, 0.020, "2026-08-05", fontsize=9, color=MUTED)

    path = OUT / "page2.jpg"
    fig.savefig(path, dpi=DPI, facecolor=BG, format="jpeg", pil_kwargs={"quality": 92})
    plt.close(fig)
    print("WROTE", path)


def main():
    page1()
    page2()
    from PIL import Image

    imgs = [Image.open(OUT / f"page{i}.jpg").convert("RGB") for i in (1, 2)]
    pdf_path = OUT / "twopager.pdf"
    imgs[0].save(pdf_path, save_all=True, append_images=imgs[1:], resolution=DPI)
    print("WROTE", pdf_path)
    for i, im in enumerate(imgs, 1):
        print(f"page{i}: {im.size[0]}x{im.size[1]} px")

    share = ROOT / "evidence" / "share"
    share.mkdir(parents=True, exist_ok=True)
    for name in ("page1.jpg", "page2.jpg", "twopager.pdf"):
        src = OUT / name
        if name.endswith(".jpg"):
            (share / f"L_closed_{name}").write_bytes(src.read_bytes())
        else:
            (share / "L_closed_two_pager.pdf").write_bytes(src.read_bytes())
    gap = 16
    w = max(imgs[0].width, imgs[1].width)
    canvas = Image.new("RGB", (w, imgs[0].height + imgs[1].height + gap), (15, 20, 25))
    canvas.paste(imgs[0], (0, 0))
    canvas.paste(imgs[1], (0, imgs[0].height + gap))
    canvas.save(share / "L_closed_two_pager_stacked.jpg", "JPEG", quality=90, optimize=True)
    canvas.save(OUT / "twopager_stacked.jpg", "JPEG", quality=90, optimize=True)
    print("mirrored to evidence/share/")


if __name__ == "__main__":
    main()
