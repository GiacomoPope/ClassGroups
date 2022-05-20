from math import gcd
from random import randint
from functools import reduce

def sign(a):
    return -1 if a < 0 else 1 

def egcd(a, b):
    sign_a, sign_b = sign(a), sign(b)
    a,b = abs(a), abs(b)
    if a == 0:
        return (b, 0, sign_b*1)
    else:
        d, v, u = egcd(b % a, a)
        return (d, sign_a*(u - (b // a) * v), sign_b*v)

def isqrt(n):
    if n == 0:
        return 0
    x = 1 << (n.bit_length() + 2 >> 1)
    while True:
        y = (x + n // x) >> 1
        if y >= x:
            return x
        x = y

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

def is_prime(n, k):
    """
    miller rabin primality test
    """
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2
    for _ in range(k):
        b = randint(2, n-1)
        x = pow(b, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def random_prime(n):
    while True:
        x = randint(2, n)
        if is_prime(x, 40):
            return x

def is_square(a, p):
    if (p & 1) == 0:
        raise Exception('Unsupported division: even modulus')
    b = p
    if a == 0:
        return True
    ls = 1
    while a != 0:
        if (a & 1) == 0:
            a >>= 1
            if ((b + 2) & 7) > 4:
                ls = -ls
        else:
            if a < b:
                a, b = b, a
                if (a & b & 3) == 3:
                    ls = -ls
            a -= b
    return ls == 1

def mod_sqrt(a, p):
    if (p & 3) == 3:
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
        raise Exception(f'Unsupported square root for this modulus: {p}')
    if pow(s,2,p) != a:
        print(a,s,pow(s,2,p))
        return False
    if (s & 1) == 0:
        s = -s % p
    return s

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