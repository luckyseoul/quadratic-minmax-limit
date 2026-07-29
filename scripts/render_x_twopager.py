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
        "Update  ·  MO 413935  ·  existence of the limit still OPEN",
        fontsize=10,
        color=MUTED,
    )
    _bar(ax, 0.905)

    _t(ax, 0.06, 0.870, "THE QUESTION  (unchanged)", fontsize=11, color=AMBER, fontweight="bold")
    _t(
        ax,
        0.50,
        0.820,
        r"$L \;=\; \lim_{n\to\infty}\, n^{-3/2}\;"
        r"\min_{a_{ij}=\pm 1}\;\max_{x_i=\pm 1}"
        r"\left|\sum_{i<j} a_{ij}\, x_i x_j\right|$",
        fontsize=13.5,
        color=WHITE,
        ha="center",
        va="center",
    )
    _t(
        ax,
        0.06,
        0.770,
        r"$\alpha_n = m_n / n^{3/2}$.  Sandwich proved; existence of $L=\lim\alpha_n$ still open.",
        fontsize=10.5,
        color=MUTED,
    )

    _card(ax, 0.05, 0.545, 0.90, 0.200)
    _t(ax, 0.08, 0.710, "PROVED  —  SANDWICH  (load-bearing)", fontsize=11, color=GREEN, fontweight="bold")
    _t(
        ax,
        0.50,
        0.645,
        r"$\frac{1}{\pi}\ \leq\ \liminf\alpha_n"
        r"\ \leq\ \limsup\alpha_n\ \leq\ \frac{1}{2}$",
        fontsize=16,
        color=WHITE,
        ha="center",
        va="center",
    )
    _t(
        ax,
        0.08,
        0.575,
        "Lower: dual-Gaussian arcsine (every Seidel matrix).\n"
        r"Upper: denseness + Paley / conference.  Does not force $\liminf=\limsup$.",
        fontsize=9.5,
        color=MUTED,
        linespacing=1.35,
    )

    _t(ax, 0.06, 0.510, "PATH TO L = 1/2  (if it exists this way)", fontsize=11, color=AMBER, fontweight="bold")
    _t(
        ax,
        0.06,
        0.470,
        r"On dense $\rho=1$ Paley family $n=p^2+1$:"
        r"  $m_n\geq\Phi-2$  $\Rightarrow$  E(1)  $\Rightarrow$  denseness  $\Rightarrow$  $L=\frac{1}{2}$.",
        fontsize=10.5,
        color=TEXT,
    )
    _t(
        ax,
        0.06,
        0.430,
        r"$\rho=1$ is proved (halfspace boolean evec $Cx=px$)."
        r"  E(1) is the blocker — reduced to residual structure on Max±.",
        fontsize=10.5,
        color=MUTED,
    )

    _card(ax, 0.05, 0.175, 0.90, 0.230)
    _t(ax, 0.08, 0.370, "WHERE WE ARE  (honest)", fontsize=11, color=AMBER, fontweight="bold")
    status_lines = [
        r"•  Sandwich + $\rho=1$ family + denseness:  shipped.",
        r"•  Bi-tight / Type-I covers:  reduced to $g_{\min}\geq L(p)$ or $\lambda_2(P\odot P)\leq 4/N$.",
        r"•  Certified $g_{\min}\geq L(p)$ at $p=5,7$  ($n=26,50$).  General $p\geq5$: OPEN.",
        r"•  Deep non-tight ND / $\Phi\geq\Phi(C)$ for $p\geq5$:  still OPEN.",
        r"•  Main Theorem $L=\frac{1}{2}$:  not claimed.  Soft-close banned.",
    ]
    y = 0.330
    for line in status_lines:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.032

    _t(
        ax,
        0.06,
        0.115,
        "MathOverflow 413935  ·  github.com/luckyseoul/quadratic-minmax-limit",
        fontsize=9,
        color=MUTED,
    )
    _t(ax, 0.94, 0.115, "1 / 2", fontsize=10, color=MUTED, ha="right")
    _t(ax, 0.06, 0.075, "Update · 2026-07-30  ·  builds on prior X cards", fontsize=9, color=MUTED)

    path = OUT / "page1.jpg"
    fig.savefig(path, dpi=DPI, facecolor=BG, format="jpeg", pil_kwargs={"quality": 92})
    plt.close(fig)
    print("WROTE", path)


def page2():
    fig, ax = _fig()

    _t(ax, 0.06, 0.955, "RESIDUAL ATTACK  ·  WHAT'S NEW", fontsize=14, color=ACCENT, fontweight="bold")
    _t(
        ax,
        0.06,
        0.922,
        "Max+ fourth moments  ·  operator T  ·  bi-tight residual  ·  still not a prize solution",
        fontsize=10,
        color=MUTED,
    )
    _bar(ax, 0.905)

    _t(ax, 0.06, 0.870, "NEW STRUCTURE  (Props 15.47–15.69, sketch)", fontsize=11, color=GREEN, fontweight="bold")

    bullets_new = [
        r"$g_{\min}=-\max|m_4|$ on 4-sets with $|\kappa|=1$  (pairing algebra).",
        r"Master identity:  $m_4=\kappa/p^2+\mathrm{Ext}/(4p)$,  $\mathrm{Ext}=Tm_4$.",
        r"Combinatorial:  $T\kappa=0$ on every $|\kappa|=1$ set  (conference $C^2$ + $K_4$).",
        r"Spectrum:  $\lambda_{\max}(T)=4p$ at $p=5,7$  $\Rightarrow$  $4pI-T$ singular;",
        r"    $m_4=m_*+h$ with $h\in\ker(4pI-T)$; min-norm $m_*$ already $\leq L$ there.",
        r"Equivalent forms:  $Q\leq16N\|B\|_F^2$  $\Leftrightarrow$  $\lambda_{\mathrm{cycle}}\leq8$"
        r"  $\Leftrightarrow$  $\lambda_2(P\odot P)\leq4/N$.",
    ]
    y = 0.830
    for b in bullets_new:
        _t(ax, 0.06, y, "▸  " + b, fontsize=9.8, color=TEXT)
        y -= 0.036

    _card(ax, 0.05, 0.455, 0.90, 0.145)
    _t(ax, 0.08, 0.565, "CERTIFIED NUMBERS  (multi-worker)", fontsize=11, color=GREEN, fontweight="bold")
    _t(
        ax,
        0.08,
        0.520,
        r"$p=5$:  $g_{\min}=-3/65 \geq L(5)=-3/50 > T(5)$."
        r"   $p=7$:  $g_{\min}=-109/2863 \geq L(7)=-5/98$.",
        fontsize=10.5,
        color=WHITE,
    )
    _t(
        ax,
        0.08,
        0.480,
        r"$p=3$ equality case for 16N / $\lambda_{\mathrm{cycle}}=8$;"
        r"  gap holds numerically $p=5,7$ (fails $p=3$).",
        fontsize=10,
        color=MUTED,
    )

    _t(ax, 0.06, 0.420, "STILL OPEN  (blockers)", fontsize=11, color=AMBER, fontweight="bold")
    open_items = [
        r"P0  Prove $g_{\min}\geq L(p)$ (or $\lambda_2(P\odot P)\leq4/N$) for all primes $p\geq5$.",
        r"P1  Deep non-tight ND:  $\Phi\geq\Phi(C)$ outside Max+-tight covers.",
        r"P2  Only if P0+P1:  Main Theorem $L=\frac{1}{2}$  (or a genuine non-existence pair).",
    ]
    y = 0.380
    for b in open_items:
        _t(ax, 0.06, y, "▸  " + b, fontsize=10.2, color=TEXT)
        y -= 0.036

    _card(ax, 0.05, 0.155, 0.90, 0.100, fc="#2a1f0a")
    _t(ax, 0.08, 0.220, "STATUS", fontsize=11, color=AMBER, fontweight="bold")
    _t(
        ax,
        0.08,
        0.175,
        r"Existence of $L$ remains OPEN.  Structure is tighter; the last mile is still proof.",
        fontsize=11,
        color=TEXT,
    )

    _bar(ax, 0.125, color="#38444d")
    _t(
        ax,
        0.06,
        0.070,
        "Handoff + code + tests:  github.com/luckyseoul/quadratic-minmax-limit\n"
        "HANDOFF.md  ·  solution.md Props 15.45–15.69  ·  evidence/  ·  not a prize solution yet",
        fontsize=9.5,
        color=MUTED,
        linespacing=1.4,
    )
    _t(ax, 0.94, 0.080, "2 / 2", fontsize=10, color=MUTED, ha="right")

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
