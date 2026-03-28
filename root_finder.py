#!/usr/bin/env python3
"""root_finder - Numerical root finding algorithms."""
import argparse, math

def bisection(f, a, b, tol=1e-12, max_iter=100):
    for i in range(max_iter):
        c = (a+b)/2; fc = f(c)
        if abs(fc) < tol or (b-a)/2 < tol: return c, i+1
        if f(a)*fc < 0: b = c
        else: a = c
    return (a+b)/2, max_iter

def newton(f, df, x0, tol=1e-12, max_iter=100):
    x = x0
    for i in range(max_iter):
        fx, dfx = f(x), df(x)
        if abs(fx) < tol: return x, i+1
        if dfx == 0: break
        x = x - fx/dfx
    return x, max_iter

def secant(f, x0, x1, tol=1e-12, max_iter=100):
    for i in range(max_iter):
        f0, f1 = f(x0), f(x1)
        if abs(f1) < tol: return x1, i+1
        if f1 == f0: break
        x2 = x1 - f1*(x1-x0)/(f1-f0)
        x0, x1 = x1, x2
    return x1, max_iter

def main():
    p = argparse.ArgumentParser(description="Root finder")
    p.add_argument("--demo", action="store_true")
    p.add_argument("-f", "--function", default="x**3-x-2", help="f(x) expression")
    p.add_argument("-a", type=float, default=1); p.add_argument("-b", type=float, default=2)
    p.add_argument("--tol", type=float, default=1e-12)
    args = p.parse_args()
    if args.demo:
        functions = [
            ("x^2 - 2 (sqrt 2)", lambda x: x**2-2, lambda x: 2*x, 1, 2),
            ("x^3 - x - 2", lambda x: x**3-x-2, lambda x: 3*x**2-1, 1, 2),
            ("cos(x) - x", lambda x: math.cos(x)-x, lambda x: -math.sin(x)-1, 0, 1),
            ("e^x - 3x", lambda x: math.exp(x)-3*x, lambda x: math.exp(x)-3, 0, 2),
        ]
        for name, f, df, a, b in functions:
            r_bi, n_bi = bisection(f, a, b)
            r_nw, n_nw = newton(f, df, (a+b)/2)
            r_sc, n_sc = secant(f, a, b)
            print(f"{name}:")
            print(f"  Bisection: x={r_bi:.12f} ({n_bi} iters)")
            print(f"  Newton:    x={r_nw:.12f} ({n_nw} iters)")
            print(f"  Secant:    x={r_sc:.12f} ({n_sc} iters)")
    else:
        f = eval(f"lambda x: {args.function}")
        r, n = bisection(f, args.a, args.b, args.tol)
        print(f"Root: {r:.15f} ({n} iterations)")

if __name__ == "__main__":
    main()
