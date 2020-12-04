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

import math

from ..plytypes import FractionFromExact, \
    FractionFromExactOrInf, fastr


class Geom1D:
    """Immutable sum of closed-open intervals of an 1D line.

    Supports the set operations ``&``, ``|``, ``^``, ``-``
    as well as ``~``..

    Basically, open intervals except that a union
    of (0, 1) and (1, 2) is (0, 2). One way to accomplish this
    theoretically is by not including the constructible numbers
    in the domain.

    So none of the exact numbers 0, 1 or 2 can ever be part of the
    infinite set of numbers contained there, only
    numbers that are greater or less, so then the middle point
    gets removed nicely.

    A simpler, equivalent definition is to use lower-bound closed,
    upper-bound open intervals a la Python ranges but the above
    is more symmetric in definition wrt the real numbers. :)

    Note that - is defined as set difference.

    The intervals can contain infinities.

            g = Geom1D([[0, 1], [2, 3]])
            g.locate(-1) --> -1
            g.locate(0) --> 0
            g.locate(0.5) --> 1
            g.locate(1.5) --> 0

            h = Geom1D([[1, 3]])

            i = g & h  # Equivalent to Geom1D([[2, 3]])
            j = g | h  # Equivalent to Geom1D([[0, 3]])
            k = g - h  # Equivalent to Geom1D([[0, 1]])
            l = g ^ h  # Equivalent to Geom1D([[0, 2]])
            m = ~h     # Equivalent to Geom1D([[-inf, 1], [3, inf]])
    """

    def __init__(self, intervals):
        intervals = sorted(intervals)

        # Check and convert
        new_intervals = []
        for interval in intervals:
            assert len(interval) == 2
            assert interval[0] <= interval[1]

            new_intervals.append([FractionFromExactOrInf(v) for v in interval])
        intervals = new_intervals

        # Normalize
        res = []
        prev = None
        for interval in intervals:
            if interval[0] == interval[1]:
                continue

            if prev is None:
                prev = interval
            else:
                if prev[1] >= interval[0]:
                    if interval[1] > prev[1]:
                        prev = (prev[0], interval[1])
                else:
                    res.append(prev)
                    prev = interval
        if prev is not None:
            res.append(prev)

        res = [tuple(r) for r in res]

        self.intervals = res

    def locate(self, point):
        """
        1 = in
        -1 = out
        0 = indeterminate, possibly on edge,
            possibly on internal, virtual edge
        """
        for interval in self.intervals:
            if interval[0] == point or interval[1] == point:
                return 0
            if interval[0] < point < interval[1]:
                return 1
            if interval[0] > point:
                break
        return -1

    def disjoint_pieces(self):
        """Return separate Geom1D pieces for disjoint intervals in self."""
        return [Geom1D([interval]) for interval in self.intervals]

    def __repr__(self):
        # fastr is not ideal but neither is the long representation...
        return "Geom1D(%s)" % (", ".join([fastr(interval)
                                          for interval in self.intervals]))

    def __fstr__(self, **args):
        return "Geom1D(%s)" % (", ".join([fastr(interval, **args)
                                          for interval in self.intervals]))

    def get_intervals(self):
        return self.intervals

    @classmethod
    def empty(self):
        """Return an empty Geom1D."""
        return Geom1D([])

    def is_empty(self):
        """Return true if this Geom1D is empty, i.e. 0-measure."""
        return len(self.intervals) == 0

    @classmethod
    def full(self):
        """Return a Geom1D that covers the full real line."""
        return Geom1D([[-math.inf, math.inf]])

    def buffer(self, amount):
        """Minkowski sum with the interval [-amount, amount].

        Amount can be negative to reduce the area.
        """
        amount = FractionFromExact(amount)
        intervals = []
        for interval in self.intervals:
            prop = [interval[0] - amount, interval[1] + amount]
            if prop[1] <= prop[0]:
                continue
            intervals.append(prop)
        return Geom1D(intervals)

    def measure1d(self):
        """Return the 1d measure (in length) of this object"""
        return sum([b - a for (a, b) in self.intervals])

    def is_bounded(self):
        """Return true if the region described by this Geom1D is bounded.
        """
        return math.isfinite(self.measure1d())

    def __or__(self, other):
        return Geom1D(self.intervals + other.intervals)

    def __and__(self, other):
        return ~((~self) | (~other))  # :)

    def __sub__(self, other):
        return self & (~other)

    def __xor__(self, other):
        return ((~self) & other) | (self & (~other))

    def __invert__(self):
        intervals = [-math.inf]
        for interval in self.intervals:
            intervals = intervals + list(interval)
        intervals = intervals + [math.inf]

        return Geom1D(zip(intervals[0::2], intervals[1::2]))

    def filter(self, f):
        """Filter all disjoint intervals in this Geom1D through the function
        and return a new Geom1D with those intervals that pass the filter.
        """
        return Geom1D([interval for interval in self.intervals
                       if f(*interval)])
