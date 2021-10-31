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

import numpy as np


def segment_segment_general_intersection(l0, l1):
    """Only when segments are in general position.
    Returns (pt, t0, t1).

    Performance is key so we don't use numpy here.
    """

    # denom = det(np.stack((l0[1] - l0[0], l1[1] - l1[0]), axis=1))
    denom = (l0[1][0] - l0[0][0]) * (l1[1][1] - l1[0][1]) - (
        l0[1][1] - l0[0][1]
    ) * (l1[1][0] - l1[0][0])

    assert denom != 0

    # Normal case, nothing to see here
    # t_num = det(np.stack([l0[0] - l1[0], l1[0] - l1[1]], axis=1))
    t_num = (l0[0][0] - l1[0][0]) * (l1[0][1] - l1[1][1]) - (
        l0[0][1] - l1[0][1]
    ) * (l1[0][0] - l1[1][0])
    # u_num = -det(np.stack([l0[0] - l0[1], l0[0] - l1[0]], axis=1))
    u_num = -(
        (l0[0][0] - l0[1][0]) * (l0[0][1] - l1[0][1])
        - (l0[0][1] - l0[1][1]) * (l0[0][0] - l1[0][0])
    )

    t0 = t_num / denom
    t1 = u_num / denom

    inter = (
        l0[0][0] + t_num * (l0[1][0] - l0[0][0]) / denom,
        l0[0][1] + t_num * (l0[1][1] - l0[0][1]) / denom,
    )

    assert inter[0] == l1[0][0] + u_num * (l1[1][0] - l1[0][0]) / denom
    assert inter[1] == l1[0][1] + u_num * (l1[1][1] - l1[0][1]) / denom

    return (inter, t0, t1)


def line_segment_point_fraction(segment, point):
    """Return the fraction the point is of the segment.

        0 = start
        ..
        1 = end
        None = not on the same line as segment
    """
    if np.all(segment[0] == segment[1]):
        return 0 if np.all(segment[0] == point) else None
    if segment[0][0] == segment[1][0]:
        # Flip x and y
        return line_segment_point_fraction(
            [(segment[0][1], segment[0][0]), (segment[1][1], segment[1][0])],
            (point[1], point[0]),
        )

    f = (point[0] - segment[0][0]) / (segment[1][0] - segment[0][0])
    if np.all(point[1] == segment[0][1] + f * (segment[1][1] - segment[0][1])):
        return f
    return None


def simple_polygon_area(poly):
    """Return the area of the simple polygon.

    CCW = positive.
    """
    area = 0

    prev = poly[-1]
    for vert in poly:
        avg_x = (vert[0] + prev[0]) / 2
        dy = -(vert[1] - prev[1])
        area += avg_x * dy
        prev = vert

    return area


def _subsign(a, b):
    """A faster way to find
    np.sign(a - b) for rationals.
    """
    if a > b:
        return 1
    if a == b:
        return 0
    return -1


def locate_point_polygon_winding_number(poly, point):
    """Locate point using winding numbers.

    Returns:
        1 if inside
        0 if on the edge
        -1 if outside

    This is currently (2020-11) one of the hottest
    functions in geom2d_py. This may change later.

    http://geomalgorithms.com/a03-_inclusion.html
    """
    winding = 0
    prev = poly[-1]
    sprev = _subsign(prev[1], point[1])
    for vert in poly:
        try:
            # Try stuff just to avoid recalculating
            # sprev
            svert = _subsign(vert[1], point[1])

            if sprev * svert == 1:
                # Edge completely above or below point
                continue

            if np.all(vert == point):
                return 0

            # Ignore horizontal unless we are on there
            if prev[1] == vert[1]:
                if (
                    prev[0] < point[0] < vert[0]
                    or prev[0] > point[0] > vert[0]
                ):
                    return 0
                continue

            # Exclude upper vertex
            if max(vert[1], prev[1]) == point[1]:
                continue

            area = simple_polygon_area([prev, vert, point])
            if area == 0:
                # On edge, exactly
                return 0
            if area > 0:
                winding += 1
            else:
                winding -= 1
        finally:
            prev = vert
            sprev = svert

    return 1 if winding != 0 else -1


def edges_iter(poly):
    """Iterate over all pairs (poly[i], poly[i + 1]), wrapping"""
    n = len(poly)
    for i in range(n):
        yield (poly[i], poly[(i + 1) % n])
