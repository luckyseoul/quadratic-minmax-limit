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
        r"MO 413935  ·  limit OPEN  ·  three load-bearing statements",
        fontsize=10,
        color=MUTED,
    )
    _bar(ax, 0.902)

    _card(ax, 0.05, 0.730, 0.90, 0.145)
    _t(ax, 0.08, 0.840, "SETUP", fontsize=11, color=AMBER, fontweight="bold")
    setup = [
        r"$m_n=\min_A\max_x |Q_A(x)|$,  $Q_A(x)=\sum_{i<j}A_{ij}x_ix_j$,  $\alpha_n=m_n/n^{3/2}$.",
        r"Question: does $\lim\alpha_n$ exist?  Identifying a value is separate.",
        r"Reviewed sandwich:  $B\leq\liminf\alpha_n\leq\limsup\alpha_n\leq 1/2$,",
        r"$B=0.3258407554\ldots>13/40$.  Convergence is not implied.",
    ]
    y = 0.800
    for line in setup:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.026

    _card(ax, 0.05, 0.380, 0.90, 0.320)
    _t(ax, 0.08, 0.665, "LEMMA A  —  paired-field lower bound  (2026-09-17)", fontsize=11, color=GREEN, fontweight="bold")
    la = [
        r"Split polynomial-phase sign covariances $C_s=C_0+s\cdot 2\kappa M/(1+q)+E_s$,",
        r"$\kappa=2/\pi$,  $\|E_s\|_F=O_L(1)$.  Retain the PSD even tail.",
        r"Pair local variances by $(M^3)_{ii}^2\leq q[(M^4)_{ii}-q^2]$.",
        r"Cap $L$ is removed by existing same-order regularization ($n\to\infty$ first).",
        r"Then $\liminf\alpha_n\geq B=[81\kappa/2+9\sqrt{\kappa(2/3+\kappa/2)}]/101$.",
        r"This is an all-orders lower bound, not a finite cutoff and not a limit.",
        r"Nine exact scalar/polynomial checks passed (NUKA).  Not a formalization.",
    ]
    y = 0.625
    for line in la:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.030

    _card(ax, 0.05, 0.175, 0.90, 0.175, fc=GREEN_CARD)
    _t(ax, 0.08, 0.310, "CHECK", fontsize=11, color=GREEN, fontweight="bold")
    chk = [
        r"•  $\kappa=2/\pi$ $\Rightarrow$ $B\approx 0.3258407554>13/40=0.325$.",
        r"•  Conference constructions still give only $\limsup\alpha_n\leq 1/2$.",
        r"•  Proof: evidence/NOTE_2026-09-17_PAIRED_POLYNOMIAL_FIELD_LOWER.md",
    ]
    y = 0.270
    for line in chk:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.028

    _t(ax, 0.06, 0.100, "github.com/luckyseoul/quadratic-minmax-limit  ·  CORE.md §4–§7", fontsize=9, color=MUTED)
    _t(ax, 0.94, 0.100, "1 / 2", fontsize=10, color=MUTED, ha="right")
    _t(ax, 0.06, 0.060, "Key lemmas only — full writeup separate", fontsize=9, color=MUTED)

    path = OUT / "lemmas_page1.jpg"
    fig.savefig(path, dpi=DPI, facecolor=BG, format="jpeg", pil_kwargs={"quality": 92})
    plt.close(fig)
    print("WROTE", path)


def page2():
    fig, ax = _fig()

    _t(ax, 0.06, 0.955, "KEY LEMMAS  (continued)", fontsize=14, color=ACCENT, fontweight="bold")
    _t(ax, 0.06, 0.920, r"two-ray  ·  residue (6.20)  ·  finite census", fontsize=10, color=MUTED)
    _bar(ax, 0.902)

    _card(ax, 0.05, 0.620, 0.90, 0.255)
    _t(ax, 0.08, 0.840, "LEMMA B  —  two-ray criterion  (CORE §7 / Prop 6.3)", fontsize=11, color=GREEN, fontweight="bold")
    lb = [
        r"Put $H(n)=m_n^{2/3}$.  If a nonnegative $\eta$ has Dini tail $E(N)\to 0$ and",
        r"    $H(2n)\leq 2H(n)+2n\eta(n)$,   $H(3n)\leq 3H(n)+3n\eta(n)$",
        r"for all large $n$, then $\alpha_n$ converges.  Proof: $h=H/n$ changes by $\eta$",
        r"along words in $\{2,3\}$; $\log 2/\log 3$ irrational $\Rightarrow$ ratio-dense semigroup.",
        r"This is sufficient, not necessary.  Neither amplification ray is proved.",
        r"Paley-skew shielding (Prop 6.6) reduces doubling to residue (6.20).",
    ]
    y = 0.800
    for line in lb:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.028

    _card(ax, 0.05, 0.340, 0.90, 0.250)
    _t(ax, 0.08, 0.555, "LEMMA C  —  residue (6.20) is live on exact minimizers", fontsize=11, color=AMBER, fontweight="bold")
    lc = [
        r"Exact $m_5=4$, $m_6=5$, $m_8=10$: residue occupies $55\%$, $91\%$, $51\%$ of pairs.",
        r"On that set, Paley-skew $R$ violates the diamond by $0.42$, $0.53$, $0.34$ in",
        r"units of $n^{3/2}$.  High-$\alpha$ random signings have sparse residue and hold.",
        r"So closing multiplier two for actual minimizers cannot use the global Paley",
        r"estimate (6.24).  Coupling of $Q_A$ with $x^TRy$ on a large set is required.",
        r"Finite only.  Does not close or refute the all-orders ray.",
    ]
    y = 0.515
    for line in lc:
        _t(ax, 0.08, y, line, fontsize=10, color=TEXT)
        y -= 0.028

    _card(ax, 0.05, 0.145, 0.90, 0.165, fc=AMBER_CARD)
    _t(ax, 0.08, 0.275, "ASSEMBLY", fontsize=11, color=AMBER, fontweight="bold")
    ass = [
        r"A + conference upper bound $\Rightarrow$ sandwich.  B is a possible close, unused.",
        r"C shows the remaining doubling gate is not an empty finite artifact.",
        r"Original limit remains OPEN.  Do not cite $L=1/2$.",
    ]
    y = 0.235
    for line in ass:
        _t(ax, 0.08, y, line, fontsize=10.2, color=TEXT)
        y -= 0.028

    _bar(ax, 0.110, color=DIM)
    _t(ax, 0.06, 0.055, "github.com/luckyseoul/quadratic-minmax-limit  ·  CORE.md / Prop 6.3, 6.6", fontsize=9.5, color=MUTED)
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
