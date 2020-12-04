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

import hypothesis as hyp
import hypothesis.strategies as hys
import numpy as np

from ....plytypes import F

from ..primitive2d import segment_segment_general_intersection


def P(a, b):
    return (F(a), F(b))


@hys.composite
def ply_fractions(draw):
    return F(draw(hys.fractions(max_denominator=101)))


def tr(transform, points):
    if transform is None:
        return points
    # return transform[:, 0:2] @ points + transform[:, 2]
    return points @ transform[:, 0:2].T + transform[:, 2]


def det22(m):
    return m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]


@hys.composite
def random_transform_matrices(draw):
    c = [draw(ply_fractions()) for _ in range(6)]

    m = np.reshape(np.array(c, dtype=object), (2, 3))
    det = det22(m)
    if det == 0:
        m[0, 0] += 1
        m[1, 1] += 1
    det = det22(m)
    hyp.assume(det != 0)

    return m


###########
#
# Test the fast general-position routine


@hyp.given(
    hys.lists(min_size=4, max_size=4,
              elements=hys.tuples(
                  ply_fractions(),
                  ply_fractions(),
              )),
    random_transform_matrices())
# Hypothesis doesn't seem to find this, hmm...
@hyp.example(
    [P(0, 0), P(0, 1), P(1, 0), P(1, 1)], np.array([[1, 0, 0], [1, 1, 0]]))
def test_segment_segment_general_intersection(pts, transform):
    def same_intersection(a, b, transform, swap, revs):
        (apt, at0, at1) = a
        (bpt, bt0, bt1) = b

        if revs[0]:
            bt0 = 1 - bt0
        if revs[1]:
            bt1 = 1 - bt1

        if swap:
            (bt0, bt1) = (bt1, bt0)

        apt = tr(transform, apt)

        if not np.all(apt == bpt):
            return (False, ('pt', a, b))

        if not (at0 == bt0):
            return (False, ('t0', at0, bt0, a, b))

        if not (at1 == bt1):
            return (False, ('t1', at1, bt1, a, b))

        return (True, '')

    def lrev(lst): return list(reversed(lst))

    op = np.array(pts, object)
    hyp.assume(det22(np.stack([op[3] - op[2], op[1] - op[0]])) != 0)

    def assrt(v, msg):
        assert v, msg

    assrt(*same_intersection(
        segment_segment_general_intersection(pts[0:2], pts[2:4]),
        segment_segment_general_intersection(pts[2:4], pts[0:2]),
        None, True, [False, False]))

    assrt(*same_intersection(
        segment_segment_general_intersection(pts[0:2], pts[2:4]),
        segment_segment_general_intersection(lrev(pts[2:4]), pts[0:2]),
        None, True, [True, False]))

    assrt(*same_intersection(
        segment_segment_general_intersection(pts[0:2], pts[2:4]),
        segment_segment_general_intersection(pts[2:4], lrev(pts[0:2])),
        None, True, [False, True]))

    assrt(*same_intersection(
        segment_segment_general_intersection(pts[0:2], pts[2:4]),
        segment_segment_general_intersection(lrev(pts[2:4]), lrev(pts[0:2])),
        None, True, [True, True]))

    tpts = pts @ transform[:, 0:2].T + transform[:, 2]

    assrt(*same_intersection(
        segment_segment_general_intersection(pts[0:2], pts[2:4]),
        segment_segment_general_intersection(tpts[0:2], tpts[2:4]),
        transform, False, [False, False]))
