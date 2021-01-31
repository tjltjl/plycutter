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

# Randomized tests for Geom2D

import hypothesis as hyp
import hypothesis.strategies as hys

from ....plytypes import F

from ..geom2d_simple_holes import Geom2D_Py_Simple_Holes as Geom2D

from .util import tr, random_transform_matrices

#########
#
# Merge squares

N = 8

squares = [[Geom2D.rectangle(x, y, x + 1, y + 1)
            for y in range(N)] for x in range(N)]
square_idxs = [(x, y) for x in range(N) for y in range(N)]


@hyp.settings(deadline=30000)
@hyp.given(
    hys.lists(elements=hys.tuples(
        hys.sampled_from(["and", "or", "sub"]),
        hys.lists(elements=hys.sampled_from(square_idxs)))),
    random_transform_matrices()
)
def test_merge_squares(ops, transform):
    cur = Geom2D.empty()
    coords = {k: 0 for k in square_idxs}
    for op, args in ops:
        if op == 'and':
            new_coords = {k: 0 for k in coords.keys()}
            for x, y in args:
                new_coords[x, y] = coords[x, y]
            coords = new_coords

        operand = Geom2D.empty()

        hyp.note(op)

        for x, y in sorted(args):
            hyp.note((x, y))
            square = squares[x][y].transformed_with(
                lambda coords: tr(transform, coords))
            operand = operand | square
            if op == 'or':
                coords[x, y] = 1
            elif op == 'sub':
                coords[x, y] = 0

        hyp.note('go')
        if op == 'and':
            cur = cur & operand
        elif op == 'or':
            cur = cur | operand
        elif op == 'sub':
            cur = cur - operand

        half = F(1, 2)
        for x, y in square_idxs:
            pt = tr(transform, (x + half, y + half))
            loc = cur.locate(pt)
            assert loc != 0
            assert (loc == 1) == (coords[x, y] == 1), coords
