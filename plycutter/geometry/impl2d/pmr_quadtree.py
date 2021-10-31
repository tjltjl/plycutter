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

from ..aabb import AABB


def round_aabb(aabb):
    return aabb


class PMRQuadTree:
    def __init__(self, aabb, intersects_aabb, items=[], threshold=4):
        """Create a new quadtree.

        Generic for different types of items.

        aabb - the main outside bounding box that must be known beforehand
        intersects_aabb - function(item, aabb) for items stored in this tree

        """
        # Round the aabb upwards
        self.aabb = round_aabb(aabb)

        self.intersects_aabb = intersects_aabb
        self.threshold = threshold

        self.root = self.Node(self, aabb)
        for item in items:
            self.add(item)

    def add(self, item):
        self.root.add(item, True)

    def find_all_similar(self, item):
        return self.find_all(lambda aabb: self.intersects_aabb(item, aabb))

    def find_all(self, intersects_aabb):
        """Find all objects that intersect certain aabbs in the tree."""
        result = set()
        self.root.find_all(intersects_aabb, result)
        return result

    class Node:
        def __init__(self, main, aabb):
            self.main = main
            self.aabb = aabb
            self.items = set()
            self.children = None

        def add(self, item, allow_split):
            if not self.main.intersects_aabb(item, self.aabb):
                return
            if self.children is None:
                if allow_split:
                    if len(self.items) + 1 >= self.main.threshold:
                        self.split()
                    allow_split = False

            if self.children is not None:
                for child in self.children:
                    child.add(item, allow_split=allow_split)
            else:
                self.items.add(item)

        def split(self):
            # Split
            aabb_mid = (self.aabb.upper + self.aabb.lower) / 2

            self.children = []

            for xr in (
                (self.aabb.lower[0], aabb_mid[0]),
                (aabb_mid[0], self.aabb.upper[0]),
            ):
                for yr in (
                    (self.aabb.lower[1], aabb_mid[1]),
                    (aabb_mid[1], self.aabb.upper[1]),
                ):
                    aabb = AABB((xr[0], yr[0]), (xr[1], yr[1]))
                    self.children.append(self.main.Node(self.main, aabb))

            # for item in self.items:
            #    for child in self.children:
            #        child.add(item, allow_split=False)

        def find_all(self, aabb_predicate, result):
            if not aabb_predicate(self.aabb):
                return
            if self.items is not None:
                result |= self.items
            if self.children is not None:
                for child in self.children:
                    child.find_all(aabb_predicate, result)


class SegmentPMRQuadTree:
    def __init__(self, aabb=None, items=[], threshold=4):
        if aabb is None:
            aabb = AABB()
            for v0, v1 in items:
                aabb.include_point(v0)
                aabb.include_point(v1)
        self.pmr = PMRQuadTree(
            aabb,
            lambda segment, aabb: aabb.intersects_segment(segment),
            items,
            threshold,
        )

    def find(self, segment):
        return self.pmr.find_all_similar(segment)
