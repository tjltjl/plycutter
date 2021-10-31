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
#
# Uses squares on a square grid to run a set of operations sequentially
# to an accumulator that starts empty.
#
# An operation looks like
#
#   ('and', [(1, 2), (3, 4)])
#
# which means that the squares at (1, 2) and (3, 4) should be 'or'ed
# together and then the accumulator should be 'and'ed with the result.

N = 8

squares = [
    [Geom2D.rectangle(x, y, x + 1, y + 1) for y in range(N)] for x in range(N)
]
square_idxs = [(x, y) for x in range(N) for y in range(N)]


@hyp.settings(deadline=30000)  # noqa: C901
@hyp.given(
    hys.lists(
        elements=hys.tuples(
            hys.sampled_from(["and", "or", "sub"]),
            hys.lists(elements=hys.sampled_from(square_idxs)),
        )
    ),
    random_transform_matrices(),
)
def test_merge_squares(ops, transform):  # noqa: C901
    # The following two are two different ways of representing
    # the same thing: accumulator as Geom2D and manual as
    # a dict of square coordinate positions.
    #
    # Later we assert that they stay the same.
    accumulator = Geom2D.empty()
    manual = {k: 0 for k in square_idxs}

    for op, args in ops:
        if op == "and":
            new_manual = {k: 0 for k in manual.keys()}
            for x, y in args:
                new_manual[x, y] = manual[x, y]
            manual = new_manual

        operand = Geom2D.empty()

        hyp.note(op)

        for x, y in sorted(args):
            hyp.note((x, y))
            square = squares[x][y].transformed_with(
                lambda coords: tr(transform, coords)
            )
            operand = operand | square
            if op == "or":
                manual[x, y] = 1
            elif op == "sub":
                manual[x, y] = 0

        hyp.note("go")
        if op == "and":
            accumulator = accumulator & operand
        elif op == "or":
            accumulator = accumulator | operand
        elif op == "sub":
            accumulator = accumulator - operand

        half = F(1, 2)
        for x, y in square_idxs:
            pt = tr(transform, (x + half, y + half))
            loc = accumulator.locate(pt)
            assert loc != 0
            assert (loc == 1) == (manual[x, y] == 1), manual
