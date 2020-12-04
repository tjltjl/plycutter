#
# plycutter - generate finger-jointed laser cutter templates from 3D objects
# Copyright (C) 2020 Tuomas J. Lukka
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

"""
Provide various types used in multiple places.

By editing this module, it is possible to, e.g., switch rational
numbers to an alternative implementation.

The reason for this class is that python very easily starts giving
out floats and we want to avoid that by converting everything
into fractions as soon as possible.
"""

import gmpy2
import math
import numpy as np
import fractions

F = Fraction = gmpy2.mpq
"""Constructor for rational numbers."""

FractionClass = gmpy2.mpq(1, 1).__class__
"""The class that rational numbers have.

Unfortunately, this is not the same as F for mpq."""

IntClass = gmpy2.mpz(1).__class__
"""The class that arbitrary-precision integers have.

Unfortunately, this is not the same as F for mpz."""

FractionFromFloat = gmpy2.f2q


def FractionFromExact(v):
    """Convert an exact value into a Fraction.

    If the value given is not exact, throw an error.
    """
    if v.__class__ == FractionClass:
        return v
    if v.__class__ == int or v.__class__ == fractions.Fraction or \
            v.__class__ == IntClass:
        return Fraction(v)
    raise Exception('Invalid input for fraction-from-exact: %s %s' %
                    (v, v.__class__))


def FractionFromExactOrInf(v):
    """Convert an exact value or an infinity into a Fraction or infinity.

    If the value given is not exact or infinity, throw an error.
    """
    if abs(v) == math.inf:
        return v
    return FractionFromExact(v)


def fastr(*objs):
    """Return a string representation where rationals are approximated.
    """
    return fstr(*objs, approx=True)


def fstr(*objs, approx=False):  # noqa (recursion)
    """Return a string representation where rationals are shown intuitively.
    """
    def recurse(o):
        if isinstance(o, np.ndarray):
            return recurse(o.tolist())
        if type(o) is tuple:
            return "(" + ", ".join([recurse(m) for m in o]) + ")"
        if type(o) is list:
            return "[" + ", ".join([recurse(m) for m in o]) + "]"
        if type(o) is set:
            return "{" + ", ".join([recurse(m) for m in o]) + "}"
        if isinstance(o, dict):
            return "{" + ", ".join(f"{recurse(k)}: {recurse(v)}"
                                   for k, v in o.items()) + "}"
        if isinstance(o, FractionClass):
            if approx:
                f = float(o)
                if math.trunc(f) == f:
                    return str(int(f))
                return str(f)
            if o.denominator == 1:
                return str(o.numerator)
            return "%s/%s" % (o.numerator, o.denominator)
        if hasattr(o, '__fstr__'):
            return o.__fstr__(approx=approx)
        return str(o)
    return "FS/%s/" % ",".join(recurse(obj) for obj in objs)


def totuple(v):
    """Convert a nested array / list / whatever nested tuples.

    Useful for creating hash keys from numpy arrays of rationals.
    """
    v = np.asarray(v)
    if len(v.shape) > 0:
        return tuple([totuple(elem) for elem in v])
    return v.tolist()
