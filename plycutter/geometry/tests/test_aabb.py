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
import hypothesis as hyp
import hypothesis.stateful
import hypothesis.strategies as hys

from ...plytypes import Fraction
F = Fraction


@hys.composite
def ply_fractions(draw):
    return F(draw(hys.fractions()))


@hys.composite
def points(draw):
    return draw(hys.tuples(ply_fractions(), ply_fractions()))


@hyp.given(
    hys.lists(elements=points(), min_size=1, max_size=10))
def test_aabb(pts):
    aabb = AABB()
    for pt in pts:
        assert not aabb.contains(pt)
    for pt in pts:
        aabb.include_point(pt)
    for pt in pts:
        assert aabb.contains(pt)
    xs = [pt[0] for pt in pts]
    ys = [pt[1] for pt in pts]
    assert not aabb.contains((xs[0], max(ys) + 1))
    assert not aabb.contains((max(xs) + 1, ys[0]))
    assert not aabb.contains((xs[0], min(ys) - 1))
    assert not aabb.contains((min(xs) - 1, ys[0]))


@hyp.given(
    hys.lists(elements=points(), min_size=1, max_size=10),
    hys.lists(elements=points(), min_size=1, max_size=10),
    hys.lists(elements=points(), min_size=1, max_size=20),
)
def test_aabb_intersection(a_pts, b_pts, test_pts):
    a = AABB()
    b = AABB()
    for pt in a_pts:
        a.include_point(pt)
    assert not b.intersects_aabb(a)
    assert not a.intersects_aabb(b)
    for pt in b_pts:
        b.include_point(pt)

    assert a.intersects_aabb(a)
    assert b.intersects_aabb(b)

    for pt in test_pts:
        if a.contains(pt) and b.contains(pt):
            assert a.intersects_aabb(b)
            assert b.intersects_aabb(a)


@hyp.given(
    hys.lists(elements=points(), min_size=1, max_size=10),
    hys.floats(min_value=0, max_value=1),
    hys.floats(min_value=0, max_value=1),
    points(),
    hys.floats(min_value=0, max_value=1e9))
def test_aabb_segment_inside(a_pts, segx, segy, to_pt, fract):
    segx = F(segx)
    segy = F(segy)
    fract = F(fract)

    aabb = AABB()
    for pt in a_pts:
        aabb.include_point(pt)

    spt = aabb.lower + (segx, segy) * (aabb.upper - aabb.lower)

    endpt = spt + fract * (spt - to_pt)

    assert aabb.intersects_segment((spt, to_pt))
    assert aabb.intersects_segment((endpt, spt))

    assert aabb.intersects_segment((endpt, to_pt))
    assert aabb.intersects_segment((to_pt, endpt))


@hyp.given(
    hys.lists(elements=points(), min_size=1, max_size=10),
    points(),
    points(),
    hys.sampled_from([0, 1]),
    hys.booleans())
def test_aabb_segment_outside(a_pts, pt1, pt2, axis, high):
    aabb = AABB()
    for pt in a_pts:
        aabb.include_point(pt)

    pt1 = [*pt1]
    pt2 = [*pt2]

    for pt in (pt1, pt2):
        if high:
            pt[axis] = max(aabb.upper[axis] + F(1, 1000000), pt[axis])
        else:
            pt[axis] = min(aabb.lower[axis] - F(1, 1000000), pt[axis])

    assert not aabb.intersects_segment((pt1, pt2))
