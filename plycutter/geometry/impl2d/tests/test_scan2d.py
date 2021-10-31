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

# from fractions import Fraction as F
from gmpy2 import mpq as F

from .. import scan2d

# Run tests using the default python fraction class to avoid extra
# dependencies. Can use the gmpy2 class equally well, the main
# code is agnostic and just requires exact numbers of some type.

######
# Vector utilities


def S(a, b):
    """Convert to segment using fractions.
    """
    return ((F(a[0]), F(a[1])), (F(b[0]), F(b[1])))


def diffv(a, b):
    return (
        a[0] - b[0],
        a[1] - b[1],
    )


def lerpv(a, b, f):
    return (
        a[0] + (b[0] - a[0]) * f,
        a[1] + (b[1] - a[1]) * f,
    )


def dotv(a, b):
    return a[0] * b[0] + a[1] * b[1]


def parallel_intersection_points(a, b):
    da = (a[1][0] - a[0][0], a[1][1] - a[0][1])
    db = (b[1][0] - b[0][0], b[1][1] - b[0][1])

    denom = scan2d._det(da, db)
    if denom != 0:
        # Not parallel
        return []

    rda = (da[1], -da[0])
    if dotv(rda, a[0]) != dotv(rda, b[0]):
        # Parallel but not collinear
        return []

    # Sort the points. Ok also when horizontal.
    spt = list(sorted([(a[0], 0), (a[1], 0), (b[0], 1), (b[1], 1),]))

    if spt[1][0] == spt[2][0]:
        # One-point overlap at middle
        return [spt[1][0]]

    if spt[0][1] == spt[1][1]:
        # No overlap (lowest two points from same segment)
        return []

    # The middle is automatically overlap
    return [spt[1][1], spt[2][1]]


def affine_transform(affine_mat, pt):
    return (
        affine_mat[0][0] * pt[0] + affine_mat[0][1] * pt[1] + affine_mat[0][2],
        affine_mat[1][0] * pt[0] + affine_mat[1][1] * pt[1] + affine_mat[1][2],
    )


######
# Hypothesis strategies


@hys.composite
def fractions(draw, *args, **kwargs):
    return F(draw(hys.fractions(*args, **kwargs)))


@hys.composite
def points(draw):
    return draw(hys.tuples(fractions(), fractions()))


@hys.composite
def proper_segments(draw):
    v0 = draw(points())
    v1 = draw(points())
    # v1 = (v1[0], v1[1] + 1)
    hyp.assume(v0 != v1)
    return (v0, v1)


@hys.composite
def proper_nonhorizontal_segments(draw):
    seg = draw(proper_segments())
    hyp.assume(seg[0][1] != seg[1][1])
    return seg


@hys.composite
def proper_affine_matrices(draw):
    """Return non-singular affine matrices."""
    m = [[draw(fractions()) for i in range(3)] for j in range(2)]

    hyp.assume(scan2d._det(m[0][0:2], m[1][0:2]) != 0)
    return m


#######
# Tests for nonparallel_intersection


@hyp.given(
    pt=points(),
    pta=points(),
    ptb=points(),
    fa=fractions(min_value=-1, max_value=0.99),
    fb=fractions(min_value=-1, max_value=0.99),
)
def test_nonparallel_intersection_point_nonparallel(pt, pta, ptb, fa, fb):
    # Generate cases from an intersection point
    # fa and fb describe the other end of the lines from pta and ptb,
    # which do intersect if
    hyp.assume(pt != pta)
    hyp.assume(pt != ptb)

    det = scan2d._det(diffv(pta, pt), diffv(ptb, pt))

    if fa <= 0 and fb <= 0 and det != 0:
        desired = pt
    else:
        desired = None

    # Determine the other ends

    pta0 = lerpv(pt, pta, fa)
    ptb0 = lerpv(pt, ptb, fb)

    # Assert all combinations

    assert (
        scan2d._nonparallel_intersection_point((pta0, pta), (ptb0, ptb))
        == desired
    )
    assert (
        scan2d._nonparallel_intersection_point((pta, pta0), (ptb0, ptb))
        == desired
    )
    assert (
        scan2d._nonparallel_intersection_point((pta0, pta), (ptb, ptb0))
        == desired
    )
    assert (
        scan2d._nonparallel_intersection_point((pta, pta0), (ptb, ptb0))
        == desired
    )


@hyp.given(
    pt=points(), pt2=points(), off=points(), f=fractions(),
)
def test_nonparallel_intersection_point_parallel(pt, pt2, off, f):
    hyp.assume(pt != pt2)
    hyp.assume(f != 0)

    # This yields two parallel segments

    seg1 = (pt, pt2)

    seg2 = (diffv(pt, off), diffv(pt2, off))
    seg2 = (seg2[0], lerpv(seg2[0], seg2[1], f))

    assert scan2d._nonparallel_intersection_point(seg1, seg2) is None


######
# The main event: tests for all_segment_intersections


def _normalize_result(intersections):
    return list(
        sorted([(pt, tuple(sorted(segs))) for pt, segs in intersections])
    )


def asi(segments):
    """Normalize inputs and outputs"""
    return _normalize_result(
        scan2d.all_segment_intersections([S(*s) for s in segments])
    )


def test_all_segment_intersections_simple():
    """Test really basic things"""

    assert asi([]) == []

    sa = S((0, 0), (1, 0))
    sb = S((0, 0), (0, 1))

    assert asi([sa, sb]) == [((F(0), F(0)), (sb, sa))]

    sa = S((-1, 0), (1, 0))
    sb = S((0, -1), (0, 1))

    assert asi([sa, sb]) == [((F(0), F(0)), (sa, sb))]


@hyp.settings(deadline=25000, suppress_health_check=[hyp.HealthCheck.too_slow])
@hyp.given(
    hys.lists(
        elements=proper_nonhorizontal_segments(), unique=True, max_size=30,
    ),
    proper_affine_matrices(),
    #    hys.lists(
    #        elements=proper_segments(),
    #        unique=True,
    #        max_size=0,  # XXX New segment endpoints
    #                     # in the middle of a collinear segment overlap
    #        # must be resolved
    #    ),
)
# , distractors):
def test_all_segment_intersections_transform(segments, transform):
    """Verify that the intersections between a bunch of segments remain the same
    under transforms and adding extra segments does not remove points"""
    distractors = []

    def transform_segment(segment):
        return (
            affine_transform(transform, segment[0]),
            affine_transform(transform, segment[1]),
        )

    tr_segments = set([transform_segment(segment) for segment in segments])

    distracted_tr_segments = list(tr_segments | set(distractors))

    intersections = asi(segments)

    distracted_tr_intersections = asi(distracted_tr_segments)

    hyp.target(
        float(sum([len(segs) for pt, segs in distracted_tr_intersections]))
    )

    manually_tr_intersections = _normalize_result(
        [
            (
                affine_transform(transform, pt),
                [transform_segment(segment) for segment in inter_segments],
            )
            for pt, inter_segments in intersections
        ]
    )

    # Filter tr_intersections
    # Fix to use walrus in python 3.8...
    tr_intersections = _normalize_result(
        [
            (pt, tr_segments & set(inter_segments))
            for pt, inter_segments in distracted_tr_intersections
            if len(tr_segments & set(inter_segments)) >= 2
        ]
    )

    # return locals()

    assert tr_intersections == manually_tr_intersections


# @hyp.given(
#         hys.lists(
#             elements=hys.tuples(
#                 proper_nonhorizontal_segments(),
#                 points())))
# def test_sweepkey_order(seg_pts):
#     # Check that the order is total
#     sweepkeys = [scan2d._SweepKey(seg, pt) for seg, pt in seg_pts]
#
#     def cmp(a, b):
#         eq = (a == b)
#         ne = (a != b)
#         lt = (a < b)
#         gt = (b < a)
#         assert eq != ne
#         if eq:
#             assert not lt
#             assert not gt
#         if ne:
#             assert lt != gt
#         if lt: return -1
#         if gt: return 1
#         return 0
#
#     # Matrix of comparison results
#     cmps = [[cmp(ka, kb) for kb in sweepkeys] for ka in sweepkeys]
#
#
#
