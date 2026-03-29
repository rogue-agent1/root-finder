#!/usr/bin/env python3
"""root_finder - Newton, bisection, secant, and Brent root finding."""
import sys, json, math

def bisection(f, a, b, tol=1e-12, maxiter=100):
    for i in range(maxiter):
        mid = (a+b)/2
        if abs(f(mid)) < tol or (b-a)/2 < tol: return mid, i+1
        if f(a)*f(mid) < 0: b = mid
        else: a = mid
    return (a+b)/2, maxiter

def newton(f, df, x0, tol=1e-12, maxiter=50):
    x = x0
    for i in range(maxiter):
        fx = f(x); dfx = df(x)
        if abs(dfx) < 1e-15: break
        x_new = x - fx/dfx
        if abs(x_new - x) < tol: return x_new, i+1
        x = x_new
    return x, maxiter

def secant(f, x0, x1, tol=1e-12, maxiter=50):
    for i in range(maxiter):
        f0, f1 = f(x0), f(x1)
        if abs(f1-f0) < 1e-15: break
        x2 = x1 - f1*(x1-x0)/(f1-f0)
        if abs(x2-x1) < tol: return x2, i+1
        x0, x1 = x1, x2
    return x1, maxiter

def brent(f, a, b, tol=1e-12, maxiter=100):
    fa, fb = f(a), f(b)
    if fa*fb > 0: return None, 0
    if abs(fa) < abs(fb): a, b = b, a; fa, fb = fb, fa
    c, fc = a, fa; d = b - a; mflag = True
    for i in range(maxiter):
        if abs(fb) < tol: return b, i+1
        if abs(fa-fc) > 1e-15 and abs(fb-fc) > 1e-15:
            s = a*fb*fc/((fa-fb)*(fa-fc))+b*fa*fc/((fb-fa)*(fb-fc))+c*fa*fb/((fc-fa)*(fc-fb))
        else:
            s = b - fb*(b-a)/(fb-fa)
        if not ((3*a+b)/4 < s < b or b < s < (3*a+b)/4): s = (a+b)/2
        fs = f(s)
        c, fc = b, fb
        if fa*fs < 0: b, fb = s, fs
        else: a, fa = s, fs
        if abs(fa) < abs(fb): a, b = b, a; fa, fb = fb, fa
    return b, maxiter

def main():
    print("Root finding demo\n")
    f = lambda x: x**3 - x - 2
    df = lambda x: 3*x**2 - 1
    for name, solver in [("Bisection", lambda: bisection(f, 1, 2)),
                          ("Newton", lambda: newton(f, df, 1.5)),
                          ("Secant", lambda: secant(f, 1, 2)),
                          ("Brent", lambda: brent(f, 1, 2))]:
        root, iters = solver()
        print(f"  {name:10s}: root={root:.12f}, f(root)={f(root):.2e}, iters={iters}")
    # Multi-root
    g = lambda x: math.sin(x) - x/3
    print(f"\nsin(x)-x/3 roots:")
    for a, b in [(-3,-2),(-.5,.5),(2,3)]:
        r, it = brent(g, a, b)
        if r is not None: print(f"  [{a},{b}]: {r:.8f}")

if __name__ == "__main__":
    main()
