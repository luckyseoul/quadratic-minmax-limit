#!/usr/bin/env python3
"""Numerical sanity check of the log-periodic model m_n = round(n^(3/2) g(ln n)).

g(t) = a + b sin(omega t) with a = 0.42, b = 0.03, omega = 0.1.  The note shows that the
integer sequence m_n is nondecreasing, satisfies m_{n+1} - m_n <= n - 2 for n >= 9 and
alpha_n = m_n / n^(3/2) in [0.39, 0.45], and that

    e_n := m_{n+1} - m_n (1 + 1/n)^(3/2) = b omega sqrt(n) cos(omega ln n) + O(1),

so e_n <= 0.003 sqrt(n) + O(1), while alpha_n = g(ln n) + O(n^(-3/2)) has no limit (the
oscillation of g has period 2 pi / omega in ln n).  This script verifies the inequalities
on 8 <= n < 2*10^6 and reports the range of alpha_n and the maximum of
e_n - b omega sqrt(n) (bounded by 1).  It cannot display the oscillation of alpha_n,
which only occurs on the scale ln n ~ 60.
Standard library only.  Run: python3 log_periodic_alpha_model_check.py
"""

import math

A, B, OMEGA = 0.42, 0.03, 0.1
NMAX = 2_000_000


def m(n):
    return round(n ** 1.5 * (A + B * math.sin(OMEGA * math.log(n))))


def main():
    prev = m(8)
    mono = True
    step = True
    amin, amax = 9.0, 0.0
    emax = -1e9
    for n in range(8, NMAX):
        nxt = m(n + 1)
        mono &= nxt >= prev
        if n >= 9:
            step &= (nxt - prev) <= n - 2
        al = prev / n ** 1.5
        amin, amax = min(amin, al), max(amax, al)
        e = nxt - prev * (1 + 1 / n) ** 1.5 - B * OMEGA * math.sqrt(n)
        emax = max(emax, e)
        prev = nxt
    print(f"8 <= n < {NMAX}: m nondecreasing: {mono}; m_(n+1)-m_n <= n-2 for n>=9: {step}")
    print(f"alpha_n in [{amin:.4f}, {amax:.4f}]")
    print(f"max over n of (e_n - b*omega*sqrt(n)) = {emax:.3f}, so e_n <= {B * OMEGA:.3f} sqrt(n) + 1 on this range")


if __name__ == "__main__":
    main()
