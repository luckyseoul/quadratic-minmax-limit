#!/usr/bin/env python3
"""
Key lemmas two-pager for X/Grok verification (no full writeup).
Output: x-cards/lemmas_page1.jpg, lemmas_page2.jpg, lemmas.pdf
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


def _bar(ax, y, color=ACCENT):
    ax.add_patch(Rectangle((0.06, y), 0.88, 0.006, facecolor=color, edgecolor="none"))


def _t(ax, x, y, s, **kw):
    defaults = dict(color=TEXT, fontfamily="sans-serif", transform=ax.transAxes)
    defaults.update(kw)
    ax.text(x, y, s, **defaults)


def page1():
    fig, ax = _fig()

    _t(ax, 0.06, 0.955, "KEY LEMMAS  (for verification)", fontsize=14, color=ACCENT, fontweight="bold")
    _t(
        ax,
        0.06,
        0.920,
        r"MO 413935  ·  $L=\lim\alpha_n=\frac{1}{2}$  ·  three load-bearing statements",
        fontsize=10,
        color=MUTED,
    )
    _bar(ax, 0.902)

    # Setup
    _card(ax, 0.05, 0.730, 0.90, 0.145)
    _t(ax, 0.08, 0.840, "SETUP", fontsize=11, color=AMBER, fontweight="bold")
    setup = [
        r"$n=p^2+1$, $p$ odd prime;  $C$ Paley conference;  $\Phi(C)=\frac{1}{2}n\sqrt{n-1}$.",
        r"Gap-2 undercutter: edge-set $F$ with $\Phi(C\oplus F)=\Phi(C)-2$.",
        r"Goal: no such $F$ descends further $\Rightarrow m_n\geq\Phi(C)-2\Rightarrow L=\frac{1}{2}$",
        r"(denseness of $\{p^2+1\}$ + universal $\limsup\alpha_n\leq\frac{1}{2}$).",
    ]
    y = 0.800
    for line in setup:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.026

    # Lemma A
    _card(ax, 0.05, 0.380, 0.90, 0.320)
    _t(ax, 0.08, 0.665, "LEMMA A  —  bi-tight covers empty  (Prop 15.167)", fontsize=11, color=GREEN, fontweight="bold")
    la = [
        r"Let $G$ be the Max$\pm$ edge Gram, $d=n/2$, mult$(\lambda_{\max})\geq d-1$, $\lambda_{\min}\geq 6$.",
        r"Majorization upper-bounds the top eigenvalue:",
        r"    $L_*=\dfrac{p^4+24p^2-1}{2(p^2-1)}$.",
        r"Algebra for primes $p\geq 5$:",
        r"    $2d-L_*=\dfrac{p^4-24p^2-1}{2(p^2-1)}>0$   (numerator $24$ at $p=5$; increasing thereafter).",
        r"Hence $\lambda_{\max}\leq L_*<2d$ $\Rightarrow$ $\lambda_{\mathrm{cycle}}<d$ $\Rightarrow$ $\lambda_{\max}(G)=n/2$ simple",
        r"$\Rightarrow$ no bi-tight size-$2p$ cover (Prop 15.55).  Residual / $16N$ not used.",
    ]
    y = 0.625
    for line in la:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.030

    # Check box
    _card(ax, 0.05, 0.175, 0.90, 0.175, fc=GREEN_CARD)
    _t(ax, 0.08, 0.310, "CHECK (algebra only)", fontsize=11, color=GREEN, fontweight="bold")
    chk = [
        r"•  At $p=5$:  $L_*=(625+600-1)/(2\cdot 24)=1224/48=25.5$,  $2d=26$,  gap $=0.5$.",
        r"•  Poly $f(x)=x^2-24x-1>0$ for $x=p^2\geq 25$.",
        r"•  Chain:  mult$+\,\lambda_{\min}$  $\Rightarrow$  $L_*$  $\Rightarrow$  $\lambda_{\mathrm{cycle}}<d$  $\Rightarrow$  bi-tight empty.",
    ]
    y = 0.270
    for line in chk:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.028

    _t(ax, 0.06, 0.100, "github.com/luckyseoul/quadratic-minmax-limit  ·  props 15.167–15.171", fontsize=9, color=MUTED)
    _t(ax, 0.94, 0.100, "1 / 2", fontsize=10, color=MUTED, ha="right")
    _t(ax, 0.06, 0.060, "Key lemmas only — full writeup separate", fontsize=9, color=MUTED)

    path = OUT / "lemmas_page1.jpg"
    fig.savefig(path, dpi=DPI, facecolor=BG, format="jpeg", pil_kwargs={"quality": 92})
    plt.close(fig)
    print("WROTE", path)


def page2():
    fig, ax = _fig()

    _t(ax, 0.06, 0.955, "KEY LEMMAS  (continued)", fontsize=14, color=ACCENT, fontweight="bold")
    _t(ax, 0.06, 0.920, r"Type I Farkas  ·  deep freeness-fail  ·  assembly", fontsize=10, color=MUTED)
    _bar(ax, 0.902)

    # Lemma B
    _card(ax, 0.05, 0.620, 0.90, 0.255)
    _t(ax, 0.08, 0.840, "LEMMA B  —  Type I freeness-fail  (Prop 15.170)", fontsize=11, color=GREEN, fontweight="bold")
    lb = [
        r"Type I, $k=3p-2$, freeness fail, gap-2 $\Rightarrow$ $s_-=-1$.  Bad case dualizes to",
        r"    $(\mathrm{Gsum}\,x)_e = 6/p-4$,   with  $0\leq x_f\leq 1$,  $\sum x=k$,  $x_e=0$.",
        r"Box–sum lower bound:  $(\mathrm{Gsum}\,x)_e \geq -12k/(p\,n)$.",
        r"Need $6/p-4 < -12k/(pn)$.  Substitute $k=3p-2$, $n=p^2+1$:",
        r"    $4p^3-6p^2-32p+18>0$ for all primes $p\geq 5$  (check $p=5$: $500-150-160+18=208>0$).",
        r"Contradiction $\Rightarrow$ freeness-fail Type I with $s_-\leq -1$ empty.",
    ]
    y = 0.800
    for line in lb:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.028

    # Lemma C
    _card(ax, 0.05, 0.340, 0.90, 0.250)
    _t(ax, 0.08, 0.555, "LEMMA C  —  deep freeness-fail  (Prop 15.171)", fontsize=11, color=GREEN, fontweight="bold")
    lc = [
        r"Deep class: $s_+=2$, freeness fail, $k\geq 3p$.",
        r"•  Free covers: weak ND (no further descent).",
        r"•  Auto-freeness for $k\leq 3p-2$; fail-equality at $k=3p-1$ empty.",
        r"•  Remaining freeness-fail dualizes to two-level equality",
        r"      $(\mathrm{Gsum}\,x)_e = 2(8-3k/p)$,  same box LB $-12k/(pn)$.",
        r"Obstructed for every prime $p\geq 5$  $\Rightarrow$ deep freeness-fail ND.",
    ]
    y = 0.515
    for line in lc:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.028

    # Assembly
    _card(ax, 0.05, 0.145, 0.90, 0.165, fc=GREEN_CARD)
    _t(ax, 0.08, 0.275, "ASSEMBLY", fontsize=11, color=GREEN, fontweight="bold")
    ass = [
        r"A + freeness / tight empty + B + C  $\Rightarrow$  $m_n\geq\Phi(C)-2$ on $n=p^2+1$.",
        r"Denseness (Prop 6.2) + $\limsup\alpha_n\leq\frac{1}{2}$  $\Rightarrow$  $L=\frac{1}{2}$.",
        r"Not claimed here: residual / $16N$ spectral package (optional; not required for $L$).",
    ]
    y = 0.235
    for line in ass:
        _t(ax, 0.08, y, line, fontsize=10.2, color=TEXT)
        y -= 0.028

    _bar(ax, 0.110, color=DIM)
    _t(ax, 0.06, 0.055, "github.com/luckyseoul/quadratic-minmax-limit  ·  props 15.167–15.171", fontsize=9.5, color=MUTED)
    _t(ax, 0.94, 0.055, "2 / 2", fontsize=10, color=MUTED, ha="right")
    _t(ax, 0.06, 0.022, "Key lemmas only — full writeup separate", fontsize=9, color=MUTED)

    path = OUT / "lemmas_page2.jpg"
    fig.savefig(path, dpi=DPI, facecolor=BG, format="jpeg", pil_kwargs={"quality": 92})
    plt.close(fig)
    print("WROTE", path)


def main():
    page1()
    page2()
    from PIL import Image

    imgs = [Image.open(OUT / f"lemmas_page{i}.jpg").convert("RGB") for i in (1, 2)]
    pdf_path = OUT / "lemmas.pdf"
    imgs[0].save(pdf_path, save_all=True, append_images=imgs[1:], resolution=DPI)
    print("WROTE", pdf_path)
    for i, im in enumerate(imgs, 1):
        print(f"lemmas_page{i}: {im.size[0]}x{im.size[1]} px")

    share = ROOT / "evidence" / "share"
    share.mkdir(parents=True, exist_ok=True)
    for name in ("lemmas_page1.jpg", "lemmas_page2.jpg", "lemmas.pdf"):
        src = OUT / name
        (share / name).write_bytes(src.read_bytes())

    gap = 16
    w = max(imgs[0].width, imgs[1].width)
    canvas = Image.new("RGB", (w, imgs[0].height + imgs[1].height + gap), (15, 20, 25))
    canvas.paste(imgs[0], (0, 0))
    canvas.paste(imgs[1], (0, imgs[0].height + gap))
    canvas.save(OUT / "lemmas_stacked.jpg", "JPEG", quality=90, optimize=True)
    canvas.save(share / "lemmas_stacked.jpg", "JPEG", quality=90, optimize=True)
    print("mirrored to evidence/share/")


if __name__ == "__main__":
    main()
