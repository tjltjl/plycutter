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
Axis-aligned bounding boxes
"""

import numpy as np
import math


def NP(v): return np.array(v, object)


class AABB:
    """Axis-aligned bounding box of points"""

    def __init__(self, lower=None, upper=None):
        self.lower = lower
        self.upper = upper
        # Use NaNs cleverly to have less special cases
        if self.lower is None:
            self.lower = [math.nan, math.nan]
        if self.upper is None:
            self.upper = [math.nan, math.nan]
        self.lower = NP(self.lower)
        self.upper = NP(self.upper)
        self._update_middle()

    def _update_middle(self):
        self.middle = (self.lower + self.upper) / 2

    def include_point(self, pt):
        """Update this AABB to include the given point.
        """
        if self.contains(pt):
            return
        for i in range(2):
            if not (pt[i] <= self.upper[i]):
                self.upper[i] = pt[i]
            if not (pt[i] >= self.lower[i]):
                self.lower[i] = pt[i]
        self._update_middle()

    def contains(self, pt):
        return np.all(pt <= self.upper) and np.all(pt >= self.lower)

    def intersects_aabb(self, other):
        return (np.all(self.upper >= other.lower)
                and np.all(other.upper >= self.lower))

    def intersects_segment(self, segment):
        """Determine whether the segment intersects this AABB.

        segment -- line segment represented as ((x0, y0), (x1, y1))
        """
        # Quick rejects
        if np.any((segment[0] < self.lower) & (segment[1] < self.lower)):
            return False
        if np.any((segment[0] > self.upper) & (segment[1] > self.upper)):
            return False

        # Separating line tests
        norm = (NP(segment[1]) - segment[0])[::-1] * (1, -1)

        mid_dot = np.dot(norm, self.middle)
        seg_dot = np.dot(norm, segment[0])

        dot_diff = mid_dot - seg_dot

        d2 = dot_diff * dot_diff

        rhs = sum(abs(norm) * (self.upper - self.middle)) ** 2

        return d2 <= rhs

    def __repr__(self):
        return f"AABB({self.lower}...{self.upper})"
