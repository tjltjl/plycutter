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


class SegmentTree1D_Slow:
    """Slow reference implementation"""

    def __init__(self, segments):
        self.segments = segments

    def stab(self, v):
        for low, high, id in self.segments:
            if low <= v and high >= v:
                yield low, high, id


class SegmentTree1D:
    """Non-incremental 1D segment tree.

    segments -- iterator with elements
        (min, max, obj)
    """

    class Node:
        def __init__(self, left, right):
            self.min = left.min
            self.closed_low = left.closed_low
            self.max = right.max
            self.closed_high = right.closed_high
            assert left.max == right.min, (left.__dict__, right.__dict__)
            assert left.closed_high == (not right.closed_low)
            self.left = left
            self.right = right
            self.cut = left.max
            self.segments = []

        def is_inside(self, pt, side=0):
            if side <= 0:
                if self.closed_high:
                    if pt > self.max:
                        return False
                else:
                    if pt >= self.max:
                        return False
            if side >= 0:
                if self.closed_low:
                    if pt < self.min:
                        return False
                else:
                    if pt <= self.min:
                        return False
            return True

        def add_segment(self, low, high, id):
            if not self.is_inside(low, side=-1):
                return
            if not self.is_inside(high, side=1):
                return

            if low <= self.min and high >= self.max:
                self.segments.append((low, high, id))
            else:
                self.left.add_segment(low, high, id)
                self.right.add_segment(low, high, id)

        def stab(self, value):
            if self.min <= value <= self.max:
                yield from iter(self.segments)
            if value < self.cut or (value == self.cut and
                                    self.left.closed_high):
                yield from self.left.stab(value)
            else:
                yield from self.right.stab(value)

        def __repr__(self):
            def c(closed): return "C" if closed else "V"
            return f"N({self.min}{c(self.closed_low)}, {self.cut}, " \
                   f"{self.max}{c(self.closed_high)}, {self.segments}, " \
                   f"{self.left}, {self.right})"

    class PointLeaf:
        def __init__(self, value):
            self.min = value
            self.max = value
            self.closed_low = True
            self.closed_high = True
            self.at = []

        def intersects_closed(self, x0, x1):
            return x0 <= self.max and x1 >= self.min

        def add_segment(self, low, high, id):
            if self.intersects_closed(low, high):
                self.at.append((low, high, id))

        def stab(self, value):
            if value == self.min:
                yield from iter(self.at)

        def __repr__(self):
            return f"PL({self.min}C, {self.max}C, {self.at})"

    class RangeLeaf:
        def __init__(self, low, high):
            self.min = low
            self.max = high
            self.closed_low = False
            self.closed_high = False
            self.at = []

        def intersects_closed(self, x0, x1):
            return x0 < self.max and x1 > self.min

        def add_segment(self, low, high, id):
            if self.intersects_closed(low, high):
                self.at.append((low, high, id))

        def stab(self, value):
            # Get both when value == self.value
            if value > self.min and value < self.max:
                yield from self.at

        def __repr__(self):
            return f"RL({self.min}V, {self.max}V, {self.at})"

    def __init__(self, segments):
        endpoints = []
        for v0, v1, _ in segments:
            assert v1 >= v0
            endpoints.extend((v0, v1))
        endpoints = list(sorted(set(endpoints)))

        if len(endpoints) == 0:
            endpoints = [0, 1]

        nodes = []
        for i in range(len(endpoints)):
            if i > 0:
                nodes.append(self.RangeLeaf(
                    endpoints[i - 1], endpoints[i]))
            nodes.append(self.PointLeaf(endpoints[i]))

        while len(nodes) > 1:
            new_nodes = []
            for i in range(0, len(nodes), 2):
                combo = nodes[i:i + 2]
                if len(combo) == 2:
                    combo = [self.Node(*combo)]
                new_nodes.extend(combo)
            nodes = new_nodes

        self.root = nodes[0]

        for segment in segments:
            self.root.add_segment(*segment)

        # print(self.root)

        assert self.root.min == min(endpoints)
        assert self.root.max == max(endpoints), (self.root.max, max(endpoints))

    def stab(self, value):
        for low, high, id in self.root.stab(value):
            assert low <= value, (value, low, high, id)
            assert high >= value, (value, low, high, id)
            yield low, high, id
