from classgroup import *
from gmpy2 import kronecker, floor, ceil, mpfr, sqrt
from math import pi, log
from sympy import factorint
from primesieve import primes

"""
Being lazy here and using imports for finding the
first n primes (for the Euler product) and for
facotring (for reducing the computed element order)
"""

DEBUG = False

def euler_product(Cl, b=10, P_bound=18):
    fps = []
    sqrt_D = sqrt(abs(Cl.D))
    P = max(2**P_bound, int(sqrt(sqrt_D)))
    sqrt_P = sqrt(mpfr(P))

    Q = sqrt_D / pi
    for p in primes(P):
        ls = kronecker(Cl.D, p)
        Qp = (1 - (ls / p))
        Q /= Qp
        # compute some random elements
        # using small primes while we 
        # loop
        if len(fps) < b and ls == 1:
            if Cl.check_prime(p):
                fps.append(Cl.lift_a(p, check=False))

    if DEBUG: print(f"Euler product: {Q}")

    B = floor(Q*(1 + 1/(2*sqrt_P)))
    C = ceil(Q*(1 - 1/(2*sqrt_P)))
    return mpz(Q), mpz(B), mpz(C), fps

def baby_steps_giant_step(g,e,B1,C1,Q1):
    # Set bounds, I increase by 0.5% each
    # side herefor good luck...
    b,c = int(1.005*B1), int(0.995*C1) 
    t = ceil((b - c) / 2)
    q = int(ceil(sqrt(t)))
    
    if DEBUG: print(f"Baby step bounds: {int(t),b,c}")

    # Baby Steps
    Id = g.parent.identity()
    x, xr = g**e, Id
    baby_steps = {}
    for r in range(q):
        baby_steps[xr] = r
        xr *= x
        if xr == Id:
            return r + 1

    """
    I cannot get Cohen's Giant steps to find a result
    Below is my own interpretation which needs to look
    at both

    Q1 ± (sq + r) ?= h

    It looks like Cohen's misses some...

    For anyone reading, alg 5.4.10 seems to ask for

    ```
    y = xr**2
    z = x**Q1
    n = Q1

    while n <= B1:
        if z in baby_steps:
            return  n - baby_steps[z]

        z_inv = z.inverse()
        
        if z_inv in baby_steps:
            return n + baby_steps[z_pos_inv]

        z *= y
        n += 2*q

    raise ValueError(f"Group order is not within the bound {B1, C1}")
    ```

    But this fails about 50% of the time.
    """

    # Compute giant steps
    # We look for s such that
    # h = Q1 ± (2qs + r)

    y = xr**2
    z_pos = x**Q1
    z_neg = z_pos.inverse()

    for s in range(q//2 + 1):
        z_pos_inv = z_pos.inverse()
        z_neg_inv = z_neg.inverse()

        if z_pos in baby_steps:
            return  Q1 + 2*s*q - baby_steps[z_pos]
        elif z_neg in baby_steps:
            return -Q1 + 2*s*q - baby_steps[z_neg]
        elif z_pos_inv in baby_steps:
            return Q1 + 2*s*q + baby_steps[z_pos_inv]
        elif z_neg_inv in baby_steps:
            return -Q1 + 2*s*q + baby_steps[z_neg_inv]

        z_pos *= y
        z_neg *= y

    raise ValueError(f"Group order is not within the bound {B1, C1}")

def reduce_element_order(g, e, n):
    Id = g.parent.identity()
    ps = factorint(n)
    for p in ps:
        for e in range(ps[p]):
            if g**(n // p) == Id:
                n = n // p
    return n

def class_number(Cl):
    Q, B, C, fps = euler_product(Cl)
    e = 1 
    B1, C1, Q1 = B, C, Q

    for g in fps:
        g = Cl.random_element(upper_bound=2**15)
        n = baby_steps_giant_step(g, e, B1, C1, Q1)
        n = reduce_element_order(g, e, n)
        e = e*n
        if e > (B - C):
            h = e*floor(B / e)
            return int(h)            
        B1 = mpz(floor(B1 / n))
        C1 = mpz(ceil(C1 / n))
    raise ValueError("Group order cannot be found, algorithm failed...")
    return None

if __name__ == '__main__':
    p = random_prime(10**15)
    Cl = ImaginaryClassGroup(-p)
    print(f"Computing: h({Cl.D})")

    h = class_number(Cl)
    print(f"h = {h}")

    score = 0
    goal = 1000
    for _ in range(goal):
        if Cl.random_element(upper_bound=2**8)**h == Cl.identity():
            score += 1
    if score != goal:
        print(f"Failed... score = {score}") 
