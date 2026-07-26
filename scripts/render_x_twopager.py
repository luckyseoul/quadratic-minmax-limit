#!/usr/bin/env python3
"""
Render a 2-page research summary as high-res JPEGs for X.

Output: x-cards/page1.jpg, page2.jpg, twopager.pdf
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


def _bar(ax, y, color=ACCENT):
    ax.add_patch(Rectangle((0.06, y), 0.88, 0.006, facecolor=color, edgecolor="none"))


def _t(ax, x, y, s, **kw):
    defaults = dict(color=TEXT, fontfamily="sans-serif", transform=ax.transAxes)
    defaults.update(kw)
    ax.text(x, y, s, **defaults)


def page1():
    fig, ax = _fig()

    _t(ax, 0.06, 0.955, "MIN–MAX  ±1  QUADRATIC  FORM", fontsize=14, color=ACCENT, fontweight="bold")
    _t(
        ax,
        0.06,
        0.922,
        "Partial results toward MO 413935  ·  existence of the limit still OPEN",
        fontsize=10,
        color=MUTED,
    )
    _bar(ax, 0.905)

    _t(ax, 0.06, 0.870, "THE QUESTION", fontsize=11, color=AMBER, fontweight="bold")
    _t(ax, 0.06, 0.835, "Does the following limit exist?", fontsize=12, color=TEXT)

    # multi-line problem (mathtext-safe)
    _t(
        ax,
        0.50,
        0.780,
        r"$L \;=\; \lim_{n\to\infty}\, n^{-3/2}\;"
        r"\min_{a_{ij}=\pm 1}\;\max_{x_i=\pm 1}"
        r"\left|\sum_{i<j} a_{ij}\, x_i x_j\right|$",
        fontsize=14,
        color=WHITE,
        ha="center",
        va="center",
    )
    _t(
        ax,
        0.06,
        0.720,
        r"Write  $\alpha_n = m_n / n^{3/2}$  with $m_n$ the min–max above.",
        fontsize=11,
        color=MUTED,
    )

    _card(ax, 0.05, 0.48, 0.90, 0.21)
    _t(ax, 0.08, 0.655, "PROVED  —  SANDWICH", fontsize=11, color=GREEN, fontweight="bold")
    _t(
        ax,
        0.50,
        0.575,
        r"$\frac{1}{\pi}\ \leq\ \liminf_{n\rightarrow\infty}\alpha_n"
        r"\ \leq\ \limsup_{n\rightarrow\infty}\alpha_n\ \leq\ \frac{1}{2}$",
        fontsize=17,
        color=WHITE,
        ha="center",
        va="center",
    )
    _t(
        ax,
        0.08,
        0.505,
        "Lower: dual-Gaussian arcsine (every Seidel matrix).\n"
        "Upper: Paley / conference matrices + denseness of Paley orders.",
        fontsize=9.5,
        color=MUTED,
        linespacing=1.35,
    )

    _t(ax, 0.06, 0.445, "EQUIVALENT FORMS", fontsize=11, color=AMBER, fontweight="bold")
    _t(
        ax,
        0.06,
        0.400,
        r"$\sum_{i<j} a_{ij} x_i x_j = \frac{1}{2}\, x^{\top} A x$"
        r"   for symmetric zero-diagonal Seidel $A$.",
        fontsize=11,
        color=TEXT,
    )
    _t(
        ax,
        0.06,
        0.355,
        r"$m_n=\min_A \Phi(A),\quad"
        r"\Phi(A)=\max_x\left|\frac{1}{2} x^{\top} A x\right|"
        r"=\frac{1}{2}\, n\, \|A\|_{\rm op}\, \rho(A)$",
        fontsize=11,
        color=TEXT,
    )
    _t(
        ax,
        0.06,
        0.310,
        r"Cut-code:  $m_n=\binom{n}{2}-2\,\rho(D_n)$"
        r"  with $D_n=\{\pm(x_i x_j)_{i<j}\}$.",
        fontsize=11,
        color=TEXT,
    )

    _card(ax, 0.05, 0.155, 0.90, 0.125, fc="#2a1f0a")
    _t(ax, 0.08, 0.245, "STATUS", fontsize=11, color=AMBER, fontweight="bold")
    _t(
        ax,
        0.08,
        0.185,
        r"Existence of $L=\lim\alpha_n$ is OPEN — neither proved nor disproved."
        "\n"
        r"The sandwich alone does not force $\liminf=\limsup$.",
        fontsize=11,
        color=TEXT,
        linespacing=1.45,
    )

    _t(
        ax,
        0.06,
        0.10,
        "MathOverflow 413935  ·  github.com/luckyseoul/quadratic-minmax-limit",
        fontsize=9,
        color=MUTED,
    )
    _t(ax, 0.94, 0.10, "1 / 2", fontsize=10, color=MUTED, ha="right")

    path = OUT / "page1.jpg"
    fig.savefig(path, dpi=DPI, facecolor=BG, format="jpeg", pil_kwargs={"quality": 92})
    plt.close(fig)
    print("WROTE", path)


def page2():
    fig, ax = _fig()

    _t(ax, 0.06, 0.955, "WHAT'S NEW  ·  WHAT'S OPEN", fontsize=14, color=ACCENT, fontweight="bold")
    _t(
        ax,
        0.06,
        0.922,
        "ρ = 1 on a dense Paley family  ·  exact small-n table  ·  blockers for existence",
        fontsize=10,
        color=MUTED,
    )
    _bar(ax, 0.905)

    _t(ax, 0.06, 0.870, "THEOREM  —  ρ = 1 FOR PALEY ORDER  p²+1", fontsize=11, color=GREEN, fontweight="bold")
    _t(
        ax,
        0.06,
        0.825,
        r"For odd prime $p$, Paley conference $C$ over $F_{p^2}$ ($n=p^2+1$)"
        r" has a halfspace boolean eigenvector:",
        fontsize=10.5,
        color=TEXT,
    )
    _t(
        ax,
        0.50,
        0.770,
        r"$C x = p\, x,\quad x\in\{\pm1\}^n"
        r"\;\;\Rightarrow\;\; \rho(C)=1,\quad"
        r"\Phi(C)=\frac{1}{2}\, n\sqrt{n-1}$",
        fontsize=13,
        color=WHITE,
        ha="center",
    )
    _t(
        ax,
        0.06,
        0.725,
        r"Along $n_k=p_k^2+1$:  $n_{k+1}/n_k\to 1$ and $\limsup\rho=1$."
        r"  Checked for $p=3,5,7$ ($n=10,26,50$).",
        fontsize=10.5,
        color=MUTED,
    )

    _t(ax, 0.06, 0.680, "EXACT  mₙ  (CERTIFIED)", fontsize=11, color=AMBER, fontweight="bold")
    headers = ["n", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    values = ["mₙ", "1", "3", "4", "4", "5", "9", "10", "12", "13"]
    xs = [0.08 + i * 0.09 for i in range(10)]
    for x, h, v in zip(xs, headers, values):
        _t(ax, x, 0.640, h, fontsize=10, color=MUTED, ha="center", fontweight="bold")
        _t(ax, x, 0.605, v, fontsize=12, color=WHITE, ha="center", fontweight="bold")
    _t(
        ax,
        0.06,
        0.560,
        r"At $n=10$ (Paley $p=3$):  $m_{10}=13 < \Phi_{\rm Paley}=15$."
        r"  Conference is not exactly optimal — E(1) must be asymptotic.",
        fontsize=10.5,
        color=TEXT,
    )

    _card(ax, 0.05, 0.28, 0.90, 0.25)
    _t(ax, 0.08, 0.495, "OPEN BLOCKERS  (FOR EXISTENCE OF L)", fontsize=11, color=AMBER, fontweight="bold")
    bullets = [
        r"E(1)  Asymptotic optimality:  $m_{n_k}=\Phi(C_k)+o(n_k^{3/2})$ along Paley.",
        r"E(2)  $\rho(C_k)\rightarrow\rho_*$ for all Paley (already $\rho=1$ on $p^2+1$).",
        r"Or: two subsequences with unequal $\lim\alpha$ (non-existence).",
        "Soft multipartite / Hadamard / Q4 / SA local-phi alone cannot finish this.",
    ]
    y = 0.450
    for b in bullets:
        _t(ax, 0.08, y, "▸  " + b, fontsize=10, color=TEXT)
        y -= 0.038

    _t(ax, 0.06, 0.240, "CONDITIONAL PICTURE", fontsize=11, color=ACCENT, fontweight="bold")
    _t(
        ax,
        0.06,
        0.195,
        r"If E(1)+E(2):  $L=\rho_*/2$.  On the $\rho=1$ family, E(1) alone"
        r" would give $\limsup\alpha=1/2$;",
        fontsize=10.5,
        color=TEXT,
    )
    _t(
        ax,
        0.06,
        0.160,
        r"matching $\liminf=1/2$ still needs more than dual-Gauss $1/\pi$."
        r"  Neither is proved.",
        fontsize=10.5,
        color=TEXT,
    )

    _bar(ax, 0.125, color="#38444d")
    _t(
        ax,
        0.06,
        0.075,
        "Handoff + code + tests:  github.com/luckyseoul/quadratic-minmax-limit\n"
        "Start at HANDOFF.md  ·  full writeup solution.md  ·  not a prize solution yet",
        fontsize=9.5,
        color=MUTED,
        linespacing=1.4,
    )
    _t(ax, 0.94, 0.085, "2 / 2", fontsize=10, color=MUTED, ha="right")

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


if __name__ == "__main__":
    main()
