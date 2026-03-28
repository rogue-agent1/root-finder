#!/usr/bin/env python3
"""Root finding: bisection, Newton-Raphson, secant, Brent."""
import math
def bisection(f,a,b,tol=1e-12,max_iter=100):
    for _ in range(max_iter):
        c=(a+b)/2
        if abs(f(c))<tol or (b-a)/2<tol: return c
        if f(c)*f(a)<0: b=c
        else: a=c
    return (a+b)/2
def newton_raphson(f,df,x0,tol=1e-12,max_iter=100):
    x=x0
    for _ in range(max_iter):
        fx=f(x);dfx=df(x)
        if abs(fx)<tol: return x
        if abs(dfx)<1e-15: break
        x-=fx/dfx
    return x
def secant(f,x0,x1,tol=1e-12,max_iter=100):
    for _ in range(max_iter):
        f0,f1=f(x0),f(x1)
        if abs(f1)<tol: return x1
        if abs(f1-f0)<1e-15: break
        x2=x1-f1*(x1-x0)/(f1-f0);x0,x1=x1,x2
    return x1
def brent(f,a,b,tol=1e-12,max_iter=100):
    fa,fb=f(a),f(b)
    if fa*fb>0: return None
    if abs(fa)<abs(fb): a,b=b,a;fa,fb=fb,fa
    c=a;fc=fa;d=b-a;mflag=True
    for _ in range(max_iter):
        if abs(fb)<tol: return b
        if fa!=fc and fb!=fc:
            s=a*fb*fc/((fa-fb)*(fa-fc))+b*fa*fc/((fb-fa)*(fb-fc))+c*fa*fb/((fc-fa)*(fc-fb))
        else: s=b-fb*(b-a)/(fb-fa)
        cond1=not((3*a+b)/4<s<b or b<s<(3*a+b)/4)
        if cond1 or abs(s-b)>=abs(d)/2: s=(a+b)/2;mflag=True
        else: mflag=False
        fs=f(s);d=b-a
        if fa*fs<0: b,fb=s,fs
        else: a,fa=s,fs
        if abs(fa)<abs(fb): a,b=b,a;fa,fb=fb,fa
        c,fc=a,fa
    return b
if __name__=="__main__":
    f=lambda x:x**3-x-2;df=lambda x:3*x*x-1
    r1=bisection(f,1,2);print(f"Bisection: {r1:.10f}")
    r2=newton_raphson(f,df,1.5);print(f"Newton: {r2:.10f}")
    r3=secant(f,1,2);print(f"Secant: {r3:.10f}")
    r4=brent(f,1,2);print(f"Brent: {r4:.10f}")
    assert all(abs(f(r))<1e-10 for r in [r1,r2,r3,r4])
    print("Root finding OK")
