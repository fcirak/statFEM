import numpy as np
from collections import deque
import baryrat

'''
- brief
Interface for rational approximation using brasil algorithm.

- details
An interface for using the BRASIL algorithm from the Baryrat library
(https://github.com/c-f-h/baryrat) for approximating the power function
x^exponent where the exponent lies in (-1,0) \cup (0,1). Given the
exponent and polynomial degree, the rational approximant is first
obtained from the brasil algorithm. The interface then converts the
barycentric form of the rational approximant to polynomial factor form.
'''

class InterfaceBrasil:
    def __init__(self, degree):
        self.m = degree
        self.exponent = 0.0
        # evaluation point for the leading polynomial coefficient
        self.eval = 0.5

    def set_exponent(self, exponent):
        self.exponent = exponent
        # relieve singularity problem for negative exponent
        self.delta = 1.0E-15 if self.exponent > 0. else 10. ** (-(5. + self.m) / 2.)

    def f(self, x):
        return x ** self.exponent

    def factors(self, polynomial, roots):
        x_minus_roots = self.eval - roots
        val = np.prod(x_minus_roots)
        coeff = polynomial(self.eval) / val
        full_factors = np.insert(roots, 0, coeff)
        return full_factors

    def rational(self):
        r = baryrat.brasil(self.f, [self.delta, 1.], [self.m, self.m])
        # obtain the numerator and its roots
        top_poly = r.numerator()
        top_roots = r.zeros()
        # workaround for baryrat's container size issue
        if (top_roots.shape[0] > self.m):
            top_roots = deque(top_roots)
            top_roots.popleft()
            top_roots = np.asarray(top_roots)
        # obtain the denominator and its roots
        bot_poly = r.denominator()
        bot_roots = r.poles()
        # full factors of numerator and denominator
        # the first element is the leading coefficient
        top_factors = self.factors(top_poly, top_roots)
        bot_factors = self.factors(bot_poly, bot_roots)
        return top_factors, bot_factors
