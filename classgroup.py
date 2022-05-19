from classgroup_helper import *

class ImaginaryClassGroup:
    """
    Here I assume that the supplied disc. is squarefree
    I can't think of a way to test for this which
    doesn't look like factoring in disguise, which 
    I don't think is worth it for cryptographically 
    large disc. Maybe worth checking around though.
    """
    def __init__(self, d):
        self.element = BinaryQuadraticForm
        if d >= -1:
            raise ValueError(f"The discriminant {d} must be smaller than -1")
        self.D = self._discriminant(d)

    def _discriminant(self, d):
        if d % 4 == 1:
            return d
        return 4*d

    def _verify_form(self, a, b, c):
        return self.D == b**2 - 4*a*c

    def _is_squarefree():
        raise NotImplementedError()

    def identity(self):
        if self.D % 4 == 0:
            return self.element(self, 1, 0)
        return self.element(self, 1, 1)

    def random_element(self):
        while True:
            a = random_prime(abs(self.D))
            D_mod_a = self.D % a
            # lazy, so i can do easy sqrt
            if (a & 3) != 3 and (a & 7) != 5:
                continue
            if is_square(D_mod_a, a):
                break
        b = mod_sqrt(D_mod_a, a)
        if self.D % 4 == 0:
            b = crt([b, 2], [(a, 1), (2, 2)], 4*a)
        sign = randint(0,1)
        if sign:
            return self.element(self, a, b)
        return self.element(self, a, -b)

    def __repr__(self):
        return f"ImaginaryClassGroup({self.D})"

    def __str__(self):
        return f"Imaginary Class Group of quadratic forms with discriminant {self.D}"

    def __call__(self, a, b, c=None):
        """
        Elements of the class group are represented as binary
        quadratic forms. Group law logic is all there in
        `BinaryQuadraticForm`. 
        """
        return self.element(self, a, b, c)


class BinaryQuadraticForm:
    def __init__(self, Cl, a, b, c=None):
        if not isinstance(Cl, ImaginaryClassGroup):
            raise TypeError("Binary Quadratic Form expects the parent to be of type `ImaginaryClassGroup`")
        self.parent = Cl

        if a < 0:
            raise ValueError(f"Binary quadratic form must be positive definite: a > 0 and D < 0. Currently: a = {a}")
        self.a = a
        self.b = b

        """
        We can either be passed `c` or compute it
        from the disc. of the Class Group. 
        """
        if c == None:
            lhs = (b**2 - self.parent.D)
            if lhs % (4*a) == 0:
                self.c = lhs // (4*a)
            else:
                raise ValueError(
                    f"The tuple {(a,b)} cannot be represented as a binary quadratic form of discriminant {self.parent.D}")
        else:
            if self.parent._verify_form(a, b, c):
                self.c = c
            else:
                raise ValueError(
                    f"The tuple {(a,b,c)} is not an binary quadratic form of discriminant {self.parent.D}")

        # This ensures that gcd(a,b,c) = 1
        self._make_primative()
        # Create a primative, positive definite binary form
        # Note: two primative elements composed is always primative
        self._reduction()

    def _to_tuple(self):
        return (self.a, self.b)

    def _make_primative(self):
        d = gcd(self.a, self.b, self.c)
        self.a //= d
        self.b //= d
        self.c //= d

    def _is_reduced(self):
        abs_b = abs(self.b)
        if abs_b < self.a and self.a < self.c:
            return True
        if abs_b == self.a and self.b >= 0:
            return True
        if self.a == self.c and self.b >= 0:
            return True
        return False

    def _reduction_euclidean_step(self):
        q, r = divmod(self.b, 2*self.a)
        if r > self.a:
            r = r - 2*self.a
            q = q + 1
        self.c -= q*(self.b + r) // 2
        self.b = r
        self._reduction_finished()

    def _reduction_finished(self):
        if self.a > self.c:
            self.b = -self.b
            self.a, self.c = self.c, self.a
            self._reduction_euclidean_step()

        elif self.a == self.c and self.b < 0:
            self.b = -self.b
        assert self._is_reduced()

    def _reduction(self):
        if -self.a < self.b and self.b <= self.a:
            self._reduction_finished()
            return
        self._reduction_euclidean_step()

    def ab(self):
        return (self.a, self.b)

    def abc(self):
        return (self.a, self.b, self.c)

    def discriminant(self):
        return self.parent.D

    def __eq__(self, other):
        if isinstance(other, BinaryQuadraticForm):
            return all([self.a == other.a, self.b == other.b, self.parent.D == other.parent.D])
        raise TypeError(
            f"{other} if of type {type(other)} and should be {type(self)}")

    def __neg__(self):
        if self.b == 0:
            return self
        return BinaryQuadraticForm(self.parent, self.a, -self.b, self.c)

    def __add__(self, other):
        if not isinstance(other, BinaryQuadraticForm):
            raise TypeError(
                f"{other} if of type {type(other)} and should be {type(self)}")
        if self.parent.D != other.parent.D:
            raise ValueError(
                "Binary forms are defined over different discriminants")

        # Step 1
        if self.a > other.a:
            a1, b1, c1 = self.a, self.b, self.c
            a2, b2, c2 = other.a, other.b, other.c
        else:
            a1, b1, c1 = self.a, self.b, self.c
            a2, b2, c2 = other.a, other.b, other.c
        s = (b1 + b2) // 2
        n = b2 - s
        # Step 2
        if a2 % a1 == 0:
            y1, d = 0, a1
        else:
            u, v, d = egcd(a2, a1)
            y1 = u
        # Step 3
        if s % d == 0:
            x2, y2, d1 = 0, -1, d
        else:
            u, v, d1 = egcd(s, d)
            x2, y2 = u, -v
        # Step 4
        v1, v2 = a1 // d1, a2 // d1
        r = (y1*y2*n - x2*c2) % v1
        a3 = v1*v2
        b3 = b2 + 2*v2*r

        return BinaryQuadraticForm(self.parent, a3, b3)

    def __iadd__(self, other):
        self = self + other
        return self

    def __sub__(self, other):
        return self + (-other)

    def __isub__(self, other):
        self = self + (-other)
        return self

    def __mul__(self, n):
        if not isinstance(n, int):
            raise TypeError(
                f"Scalar multiplication must be performed using an integer.")

        Q = self
        R = self.parent.identity()
        # Deal with negative scalar multiplication
        if n < 0:
            n = -n
            Q = -Q
        while n > 0:
            if n % 2 == 1:
                R = R + Q
            Q = Q + Q
            n = n // 2
        return R

    def __rmul__(self, other):
        return self * other

    def __hash__(self):
        return hash(self._to_tuple())

    def __repr__(self):
        return f"BinaryQuadraticForm{self.a, self.b, self.c}"

    def __str__(self):
        return f"Binary Quadratic form {self.a, self.b} with discriminant {self.parent.D}"
