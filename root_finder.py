#!/usr/bin/env python3
"""root_finder - Bisection, Newton, secant, Brent root finding."""
import sys, argparse, math

def bisection(f, a, b, tol=1e-10, max_iter=100):
    if f(a) * f(b) > 0: raise ValueError("f(a) and f(b) must have different signs")
    for i in range(max_iter):
        c = (a + b) / 2
        if abs(f(c)) < tol or (b - a) / 2 < tol: return c, i+1
        if f(a) * f(c) < 0: b = c
        else: a = c
    return (a + b) / 2, max_iter

def newton(f, df, x0, tol=1e-10, max_iter=100):
    x = x0
    for i in range(max_iter):
        fx = f(x); dfx = df(x)
        if abs(dfx) < 1e-15: raise ValueError("Zero derivative")
        x_new = x - fx / dfx
        if abs(x_new - x) < tol: return x_new, i+1
        x = x_new
    return x, max_iter

def secant(f, x0, x1, tol=1e-10, max_iter=100):
    for i in range(max_iter):
        f0, f1 = f(x0), f(x1)
        if abs(f1 - f0) < 1e-15: break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x2 - x1) < tol: return x2, i+1
        x0, x1 = x1, x2
    return x1, max_iter

def brent(f, a, b, tol=1e-10, max_iter=100):
    fa, fb = f(a), f(b)
    if fa * fb > 0: raise ValueError("f(a) and f(b) must have different signs")
    if abs(fa) < abs(fb): a, b = b, a; fa, fb = fb, fa
    c, fc = a, fa; d = b - a; mflag = True
    for i in range(max_iter):
        if abs(fb) < tol: return b, i+1
        if abs(fa - fc) > 1e-15 and abs(fb - fc) > 1e-15:
            s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb))
        else:
            s = b - fb*(b-a)/(fb-fa)
        conds = [not ((3*a+b)/4 < s < b or b < s < (3*a+b)/4),
                 mflag and abs(s-b) >= abs(b-c)/2,
                 not mflag and abs(s-b) >= abs(c-d)/2]
        if any(conds): s = (a+b)/2; mflag = True
        else: mflag = False
        fs = f(s); d = c; c = b; fc = fb
        if fa * fs < 0: b = s; fb = fs
        else: a = s; fa = fs
        if abs(fa) < abs(fb): a, b = b, a; fa, fb = fb, fa
    return b, max_iter

def main():
    p = argparse.ArgumentParser(description="Root finding algorithms")
    p.add_argument("--demo", action="store_true")
    args = p.parse_args()
    if args.demo:
        f = lambda x: x**3 - x - 2
        df = lambda x: 3*x**2 - 1
        print("Finding root of x³ - x - 2 = 0")
        for name, solver in [("Bisection", lambda: bisection(f, 1, 2)),
                              ("Newton", lambda: newton(f, df, 1.5)),
                              ("Secant", lambda: secant(f, 1, 2)),
                              ("Brent", lambda: brent(f, 1, 2))]:
            root, iters = solver()
            print(f"  {name:12s}: x={root:.10f} f(x)={f(root):.2e} iters={iters}")

        print("\nFinding sqrt(2) via x² - 2 = 0")
        g = lambda x: x*x - 2
        root, _ = brent(g, 1, 2)
        print(f"  sqrt(2) = {root:.15f} (exact: {math.sqrt(2):.15f})")
    else: p.print_help()

if __name__ == "__main__":
    main()
