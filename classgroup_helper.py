from gmpy2 import gcd, gcdext, isqrt, is_prime, mpz, kronecker
from random import randint
from functools import reduce

def egcd(a, b):
    return gcdext(a,b)

def part_eucl(a, b, L):
    """
    Sub algorithm from Cohen p.248
    """
    # [Initialize]
    v, v2, v3 = 0, 1, b
    d, z = a, 0
    # [Euclidean Step]
    while abs(v3) > L:
        q, t3 = divmod(d, v3)
        # ensure 0 <= t3 < |v3|
        if t3 < 0:
            q  += 1
            t3 -= v3
        t2 = v - q*v2
        v, d = v2, v3
        v2, v3 = t2, t3
        z = z + 1
    # [Finished]      
    if z % 2 == 1:
        v2, v3 = -v2, -v3
    return z, d, v, v2, v3

def random_prime(n):
    while True:
        x = randint(2, n)
        if is_prime(x):
            return x

def is_square(a, p):
    return kronecker(a, p) == 1

def mod_sqrt(a, p):
    if p == 2:
        return a % p
    elif (p & 3) == 3:
        s = pow(a, (p+1) // 4, p)
    elif (p & 7) == 5:
        # Atkin's formulas:
        #   b <- (2*a)^((m-5)/8)
        #   c <- 2*a*b^2
        #   return a*b*(c - 1)
        b = pow(2*a, ((p - 5) // 8), p)
        c = 2*a*b**2
        s = a*b*(c - 1) % p
    else:
        s = tonelli_shanks(a, p)
    if pow(s,2,p) != a:
        return False
    if (s & 1) == 0:
        s = -s % p
    return s

def tonelli_shanks(a, p):
    """
    Alg. 1.5.1 Cohen (page 33)
    """
    def reduce_prime(p):
        q = p - 1
        e = 0
        while q % 2 == 0:
            q //= 2
            e += 1
        return q, e

    def find_generator(q, p):
        while True:
            n = randint(0, p)
            if kronecker(n, p) == -1:
                return pow(n, q, p)

    def find_exponent(b):
        bm, m = b, 1
        while True:
            bm = pow(bm, 2, p)
            if bm == 1:
                return m
            m += 1

    # Initialise
    q, r = reduce_prime(p)
    y = find_generator(q, p)
    x = pow(a, (q-1) // 2, p)
    b = (a*x**2) % p
    x = (a*x) % p

    #Find and reduce exp
    while b != 1:
        m = find_exponent(b) 
        if m == r:
            raise ValueError(f'a is not a quadratic residue modulo p: {a, p}')

        t_exp = pow(2, r - m - 1)
        t = pow(y, t_exp, p)
        y = pow(t, 2, p)
        r = m
        x = (x*t) % p
        b = (b*y) % p

    return x

def crt(xs, ns_fac, n):
    x = 0
    ns = [p ** e for p, e in ns_fac]
    common = reduce(gcd, ns)
    ns = [n // common for n in ns]

    for xi, ni in zip(xs, ns):
        yi = n // ni
        zi = pow(yi, -1, ni)
        x += xi * yi * zi
    return x % n