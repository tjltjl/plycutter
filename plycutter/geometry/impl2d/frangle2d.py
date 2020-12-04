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

from ..types2d import Coord, Point, Vector

"""Fractional angles: unit-circle analogues of radians

Useful for for example sorting.
"""


MAX_FRANGLE = 8
"""Frangle corresponding to full circle (square)"""

FRANGLE_90 = 2
"""Frangle corresponding to 90 degrees.

Unlike most frangles, this is conveniently invariant:
adding this to a frangle will rotate it by exactly 90
degrees."""

FRANGLE_180 = 4
"""Frangle corresponding to 180 degrees.

Unlike most frangles, this is conveniently invariant:
adding this to a frangle will rotate it by exactly 90
degrees."""


def vector_frangle(vec: Point) -> Coord:
    """An fractional angle for order key for line; sorting by these
    will yield a clockwise order starting from vertical line.

    Kind of like an angle but rational, goes from 0 to 8 along
    the (-1, -1) .. (1, 1) square

    Properties of frangles:

        fa - fb == 8 --> 360 degrees
        fa - fb == 4 --> 180 degrees
        fa - fb == 2 --> 90 degrees

    For 45 degrees, the formula would not be constant
    """
    assert vec[0] != 0 or vec[1] != 0

    if abs(vec[0]) > abs(vec[1]):
        order = abs(vec[1] / vec[0])
    else:
        order = 2 - abs(vec[0] / vec[1])

    if vec[0] >= 0:
        if vec[1] <= 0:
            return order
        else:
            return 8 - order
    else:
        if vec[1] <= 0:
            return 4 - order
        else:
            return 4 + order


def frangle_unit_square_vector(frangle: Coord) -> Vector:
    """Get a vector at frangle on the unit square.
    """
    frangle = frangle % 8

    frhalf = frangle % 4

    frquad = abs(2 - frhalf)

    if frquad < 1:
        v = (frquad, -1)
    else:
        v = (1, -(2 - frquad))

    if frhalf > 2:
        v = (-v[0], v[1])

    if frangle >= 4:
        v = (-v[0], -v[1])

    return v
