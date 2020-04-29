from random import random
from math import *


def powVar(min, max, n):
    """Returns a random number with a power law distribution of power n,      
    from values between min and max."""

    range = pow(max, 1+n) - pow(min, 1+n)
    return pow(random()*range + pow(min, 1+n), 1.0/(1+n))



def randM(min, max, a):
    """Generates a random value for the mass of a star
    where min and max define the range for the value of
    the mass and a is the exponent associated with the
    probability distribution for masses between 1 and 150
    solar masses. min can be no smaller than 0.08 and max
    can be no more than 150. a must be within the range
    2.35-3.2."""

    if max < min:
        raise ValueError("minM must be less than or equal to maxM")

    if min < .08:
        raise ValueError("minM must be no less than .08.")

    if max > 150:
        raise ValueError("maxM must be no more than 150.")

    if a < 2.35 or a > 3.2:
        raise ValueError("aIMF must be between 2.35 and 3.2")

    coef = pow(.5, -2.2)/pow(.5, -1.3)
    p1 = -1.3
    p2 = -2.2

    c = 1.87
    u = 1.0
    v = 1.0

    def f(x):
        """The probability distribution function"""
        if x >= .08 and x < .5:
            return coef*pow(x, p1)
        elif x >= .5  and x < 1.0:
            return pow(x, p2)
        else:
            return pow(x, -a)

    def g(x):
        """The comparison function"""
        return pow(x, p1)

    while c*u > f(v)/g(v):
        # reject v and pick a new one
        u = random()
        v = powVar(min, max, -1.3)

    return v



def randR(min, max):
    """Returns a rondom number for the radius of separation
    of a binary. min cannot be zero"""

    if min <= 0:
        raise ValueError("min cannot be 0 or negative")

    if max < min:
        raise ValueError("min must be less than or equal to max")

    y = random()*(log(max) - log(min))
    return exp(y)*min





