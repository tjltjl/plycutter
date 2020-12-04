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
from ..geom2d import Geom2D
from ...plytypes import F


def geom(n, offset=(0, 0)):
    pts = []
    for i in range(n):
        pts.append((F(i) + offset[0], F(i * i) + offset[1]))
    return Geom2D.polygon(pts, True)


@pytest.mark.parametrize('n', [12, 24, 48, 96])
@pytest.mark.benchmark()
def test_intersection(n, benchmark):
    a = geom(n)
    b = geom(n, offset=(F(1, 10), F(1, 10)))
    benchmark(lambda: a & b)


@pytest.mark.parametrize('n', [12, 24, 48, 96])
@pytest.mark.benchmark()
def test_union(n, benchmark):
    a = geom(n)
    b = geom(n, offset=(F(1, 10), F(1, 10)))
    benchmark(lambda: a | b)


@pytest.mark.parametrize('n', [12, 24, 48, 96])
@pytest.mark.benchmark()
def test_dilate(n, benchmark):
    a = geom(n)
    benchmark(lambda: a.buffer(F(1, 10)))


@pytest.mark.parametrize('n', [12, 24, 48, 96])
@pytest.mark.benchmark()
def test_erode(n, benchmark):
    a = geom(n)
    benchmark(lambda: a.buffer(F(-1, 10)))
