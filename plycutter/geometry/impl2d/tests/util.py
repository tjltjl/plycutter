#
# plycutter - generate finger-jointed laser cutter templates from 3D objects
# Copyright (C) 2021 Tuomas J. Lukka
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

# Utilities for tests

import hypothesis as hyp
import hypothesis.strategies as hys
import numpy as np

from ....plytypes import F


def tr(transform, points):
    """Transform the 2D points using the (2, 3) matrix"""
    if transform is None:
        return points
    # return transform[:, 0:2] @ points + transform[:, 2]
    return points @ transform[:, 0:2].T + transform[:, 2]


@hys.composite
def ply_fractions(draw):
    """Draw fractions suitable for use with geometry.

    Keeping the max_denominator smallish to keep runtime low."""
    return F(draw(hys.fractions(max_denominator=101)))


def det22(m):
    """Determinant of 2x2 matrix"""
    return m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]


@hys.composite
def random_transform_matrices(draw):
    """Random non-singular transform matrices"""
    c = [draw(ply_fractions()) for _ in range(6)]

    m = np.reshape(np.array(c, dtype=object), (2, 3))
    det = det22(m)
    if det == 0:
        m[0, 0] += 1
        m[1, 1] += 1
    det = det22(m)
    hyp.assume(det != 0)

    return m
