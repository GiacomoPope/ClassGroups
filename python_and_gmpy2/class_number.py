from classgroup_gmpy2 import *
from gmpy2 import kronecker, floor, ceil, mpfr, sqrt
from math import pi
from sympy import factorint

def primes_to_n(n):
    # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    # Returns  a list of primes < n
    sieve = [True] * n
    for i in range(3, int(n ** 0.5) + 1, 2):
        if sieve[i]:
            sieve[i * i::2 * i] = [False] * ((n - i * i - 1) // (2 * i) + 1)
    return [2] + [i for i in range(3, n, 2) if sieve[i]]

def euler_product(Cl, b=10):
    fps = []
    sqrt_D = sqrt(mpfr(abs(Cl.D)))
    P = max(2**18, sqrt(sqrt_D))
    sqrt_P = sqrt(mpfr(P))

    Q = sqrt_D / mpfr(pi)
    for p in primes_to_n(P):
        ls = kronecker(Cl.D, p)
        Qp = (1 - ls / mpfr(p))
        Q /= Qp
        # compute some random elements
        # using small primes while we 
        # loop
        if len(fps) < b and ls == 1:
            if Cl.check_prime(p):
                fps.append(Cl.lift_a(p, check=False))

    B = floor(Q*(1 + 1/(2*sqrt_P)))
    C = ceil(Q*(1 - 1/(2*sqrt_P)))
    return mpz(Q), mpz(B), mpz(C), fps

def baby_steps_giant_step(g,e,B1,C1,Q1):
    baby_steps = {}
    q = int(ceil(sqrt((B1 - C1)/ 2)))
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
            return r
    # Organise...
    y = x1 + xr
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

def class_number(Cl):
    Q, B, C, fps = euler_product(Cl)
    baby_steps = {}
    e = 1 
    B1 = B
    C1 = C
    Q1 = Q

    for g in fps:
        n = baby_steps_giant_step(g, e, B1, C1, Q1)
        n = reduce_element_order(g, e, n)
        e = e*n
        if e > (B - C):
            h = e*floor(B / e)
            return int(h)            
        B1 = floor(B1 / n)
        C1 = ceil(C1 / n)
    raise ValueError("Group order cannot be found, algorithm failed...")
    return None

p = random_prime(10**10)
Cl = ImaginaryClassGroup(-p)
h = class_number(Cl)
print(f"h(-{p}) = {h}")

score = 0
for _ in range(50):
    if h*Cl.random_element(upper_bound=2**16) == Cl.identity():
        score += 1
print(f"score: {score}")

            
