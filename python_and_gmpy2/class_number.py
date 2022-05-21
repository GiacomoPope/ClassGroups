from classgroup_gmpy2 import *
from gmpy2 import kronecker, floor, ceil, mpfr, sqrt
from math import pi, log
from sympy import factorint
from primesieve import primes

def euler_product(Cl, b=5, p_bound=18):
    fps = []
    sqrt_D = sqrt(abs(Cl.D))
    P = max(2**p_bound, int(sqrt(sqrt_D)))
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

    print(f"Euler product: {Q}")
    B = floor(Q*(1 + 1/(2*sqrt_P)))
    C = ceil(Q*(1 - 1/(2*sqrt_P)))
    return mpz(Q), mpz(B), mpz(C), fps

def baby_steps_giant_step(g,e,B1,C1,Q1):
    baby_steps = {}
    q = int(ceil(sqrt((B1 - C1)/2)))
    print(f"Baby step bounds: {q, B1, C1}")
    # Baby Steps
    Id = g.parent.identity()
    x1 = e*g
    xr = x1
    baby_steps[Id] = 0
    baby_steps[x1] = 1
    for r in range(2, q):
        xr = x1 + xr
        baby_steps[x1] = r
        if xr == Id:
            print("tada!")
            return r
    # Prepare
    if q > 2:
        y = x1 + xr
    else:
        y = q*x1
    y = 2*y
    z = Q1*x1
    n = Q1
    # Giant Steps
    while n <= B1:
        if z in baby_steps:
            return n - baby_steps[z]
        z_inv = -z
        if z_inv in baby_steps:
            return n + baby_steps[z_inv]
        z = z + y
        n += 2*q
    raise ValueError("Group order is larger than the bound {B}")

def reduce_element_order(g, e, n):
    Id = g.parent.identity()
    ps = factorint(n)
    for p in ps:
        for e in range(ps[p]):
            if (n // p)*g == Id:
                n = n // p
    return n

def class_number(Cl, p_bound=18):
    Q, B, C, fps = euler_product(Cl, p_bound=p_bound)
    e = 1 
    B1, C1, Q1 = B, C, Q

    for g in fps:
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

p = random_prime(10**8)
Cl = ImaginaryClassGroup(-p)
h = class_number(Cl, p_bound=22)
print(f"h(-{p}) = {h}")

score = 0

print(f"Testing 100 random elements to see if h*g = Id...")
for _ in range(100):
    if h*Cl.random_element(upper_bound=2**16) == Cl.identity():
        score += 1
print(f"score: {score}")

            
