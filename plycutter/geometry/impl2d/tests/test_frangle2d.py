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

import pytest
import hypothesis as hyp
import hypothesis.stateful
import hypothesis.strategies as hys
import numpy as np

from ....plytypes import F, fstr
from ..frangle2d import vector_frangle, frangle_unit_square_vector


@hys.composite
def ply_fractions(draw):
    return F(draw(hys.fractions(max_denominator=101)))


@hys.composite
def points(draw):
    return draw(hys.tuples(ply_fractions(), ply_fractions()))


def P(a, b):
    return (F(a), F(b))


@pytest.mark.parametrize("fx,fy,o", [
    (4, 0, 0),
    (4, -2, 1/2),
    (4, -4, 1),
    (2, -4, 3/2),
    (0, -4, 2),
    (-2, -4, 5/2),
    (-4, -4, 3),
    (-4, -2, 7/2),
    (-4, 0, 4),
    (-4, 2, 9/2),
    (-4, 4, 5),
    (-2, 4, 11/2),
    (0, 4, 6),
    (2, 4, 13/2),
    (4, 4, 7),
    (4, 2, 15/2),
])
def test_vector_frangle(fx, fy, o):
    assert vector_frangle((F(fx, 16), F(fy, 16), None)) == o
    assert vector_frangle(frangle_unit_square_vector(
        vector_frangle((F(fx, 16), F(fy, 16), None)))) == o


@hyp.given(points())
@hyp.example(P(1, 0))
def test_vector_frangle_2(pt):
    pt = np.array(pt, object)
    hyp.assume(sum(pt ** 2) > 0)

    frangle = vector_frangle(pt)
    usq = frangle_unit_square_vector(frangle)

    frangle2 = vector_frangle(usq)

    hyp.note(fstr(pt, frangle, usq, frangle2))

    assert frangle == frangle2

    assert usq[0] * pt[1] == usq[1] * pt[0]
    assert (usq[0] == 0) == (pt[0] == 0)
    assert (usq[1] == 0) == (pt[1] == 0)
