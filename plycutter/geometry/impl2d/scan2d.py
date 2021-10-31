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

"""2D line sweep scanning algorithms
"""
import math
from sortedcontainers import SortedList, SortedDict

# Setting this to True enables some expensive
# checks used in the tests
_extra_checks = False


def all_segment_intersections(segments):
    """Find all intersections between the segments.

    Returns iterator over entries (pt, (segment, ...))

    The segments must be unique.

    If two segments overlap on a line, only the ends of the
    overlap create intersection points.

    If given segments that use the fractions.Fraction class
    or the gmpy.mpq class or some other exact rational arithmetics.
    Due to the way Python works, make sure that all inputs
    are of the rational class, even if some coordinates are
    representable as integers such as 0, 1.

    The segments can be floating-point in which case the results
    will be inexact but should still work (though this case
    has not been tested).

    Internals:

    Using the Bentley-Ottmann algorithm,
    adapted from Berg et al: Computational Geometry, 2nd edition

    However, instead of the extra complexity of handling
    horizontal segments, we simply find an explicit transformation

        y <- y + alpha x

    which leaves the points in general position.
    """

    segments = list(segments)
    alpha = _find_general_y_shear(segments)

    tr_segments = [
        tuple([(pt[0], pt[1] + alpha * pt[0]) for pt in segment])
        for segment in segments
    ]
    origs = {k: v for k, v in zip(tr_segments, segments)}

    for pt, segs in _all_segment_intersections_no_horizontal(tr_segments):
        npt = (pt[0], pt[1] - alpha * pt[0])
        assert type(npt[0]) != float
        assert type(npt[1]) != float
        yield (npt, [origs[seg] for seg in segs])


def _find_general_y_shear(segments):
    """Find a transformation so no segments are horizontal.

    Could do something nicer to get less decimals
    but for now, just take alpha to be zero if no horizontals
    or half of the smallest non-horizontal direction if there are any.

    Just to avoid sign errors, we take the minimum over slopes
    in both directions.
    """
    had_horizontal = False
    max_ratio = None

    for segment in segments:
        d = (segment[1][0] - segment[0][0], segment[1][1] - segment[0][1])
        if d[0] != 0:
            if d[1] == 0:
                had_horizontal = True
            else:
                ratio = abs(d[1] / d[0])
                if max_ratio is None or ratio < max_ratio:
                    max_ratio = ratio
    if not had_horizontal:
        return 0
    if max_ratio is None:
        return 1
    return max_ratio / 2


def _at_y_no_horizontal(segment, y):
    if y == segment[0][1]:
        return segment[0][0]
    if y == segment[1][1]:
        return segment[1][0]
    # TODO cache
    dy = segment[1][1] - segment[0][1]
    assert dy != 0
    x0 = segment[0][0]
    dx = segment[1][0] - segment[0][0]
    dy_nw = y - segment[0][1]
    return x0 + dx * dy_nw / dy


class _SweepKey:
    """A key for the y line sweep.

    The key takes into account where the segment is inserted
    and compares the segmments at and after that point.

    This avoids the complexity of inserting things into a binary
    tree with different orders and breaking the encapsulation
    between tree balancing and the sweep logic.

    The keys are not totally ordered but the subset that is
    in the active data structure at one time is.
    """

    def __init__(self, segment, pt):
        # Horizontal is not allowed
        assert segment[0][1] != segment[1][1]
        self.segment = segment
        self.pt = pt

    def __repr__(self):
        return f"_SweepKey({self.segment}, {self.pt})"

    def __hash__(self):
        return hash((self.segment, self.pt))

    def __eq__(self, other):
        return (self.pt == other.pt) and (self.segment == other.segment)

    def __lt__(self, other):
        # We compare at the maximum event y
        if other.pt[1] > self.pt[1]:
            cmp_y = other.pt[1]
            other_x = other.pt[0]
            self_x = self.at_y(cmp_y)
        else:
            cmp_y = self.pt[1]
            other_x = other.at_y(cmp_y)
            self_x = self.pt[0]

        if self_x < other_x:
            return True
        if self_x > other_x:
            return False

        # Now we know that self_x == other_x, i.e., the lines
        # meet at the relevant y.

        # If the points are from the same Y, must compare slightly after;
        # if from different Ys, must compare slightly before
        # to get results consistent with others
        if other.pt[1] == self.pt[1]:
            alt_cmp_y = cmp_y + 1
        else:
            alt_cmp_y = cmp_y - 1
        self_next_x = self.at_y(alt_cmp_y)
        other_next_x = other.at_y(alt_cmp_y)

        if self_next_x < other_next_x:
            return True
        if self_next_x > other_next_x:
            return False

        if self.segment < other.segment:
            return True

        return False

    def at_y(self, y):
        return _at_y_no_horizontal(self.segment, y)


def _all_segment_intersections_no_horizontal(segments):  # noqa
    # Must be unique
    assert len(set(segments)) == len(segments)
    segments = list(segments)

    # Must not be degenerate
    for segment in segments:
        assert segment[0] != segment[1]

    # Use the convention from the book: sweep on Y axis
    def event_key(pt):
        return (pt[1], pt[0])

    # From point to list of segments
    event_queue = SortedDict(event_key)

    def add_event(pt, segment_key=None):
        if pt not in event_queue:
            event_queue[pt] = []
        if segment_key is not None:
            event_queue[pt].append(segment_key)

    for i, segment in enumerate(segments):
        if event_key(segment[0]) < event_key(segment[1]):
            add_event(segment[0], _SweepKey(segment, segment[0]))
            add_event(segment[1], None)
        else:
            add_event(segment[0], None)
            add_event(segment[1], _SweepKey(segment, segment[1]))

    active = SortedList()

    y = -math.inf

    while len(event_queue) > 0:
        v = event_queue.popitem(0)
        pt, segstarts = v

        # Can't be > since while there are no horizontal segments,
        # there can still be points in horizontal relation to one another
        assert pt[1] >= y
        y = pt[1]

        # Find all segments within the event point

        fake_segment = ((pt[0], pt[1]), (pt[0], pt[1] + 1))
        fake_key = _SweepKey(fake_segment, pt)

        touches = []

        # The next lower / higher keys, respectively, to enter new events for
        neighbours = []

        if _extra_checks:
            _assert_fully_sorted(list(active), y)
        # Iterate on both sides
        for it in (
            active.irange(
                None, fake_key, inclusive=(True, True), reverse=True
            ),
            active.irange(fake_key, None, inclusive=(False, True)),
        ):
            neighbour = None
            for sweep_key in it:
                if sweep_key.at_y(y) != pt[0]:
                    neighbour = sweep_key
                    break
                touches.append(sweep_key)
            neighbours.append(neighbour)

        # Remove the old sweep keys
        for touch in touches:
            active.remove(touch)

        segments_at_pt = [
            sweep_key.segment for sweep_key in touches + segstarts
        ]
        if len(segments_at_pt) > 1:
            yield (pt, tuple(segments_at_pt))

        # Create new _SweepKeys, automatically sorts
        # according to order after point
        sweep_keys = []
        for segment in segments_at_pt:
            # Is this segment still relevant?
            if max(segment[0][1], segment[1][1]) <= pt[1]:
                continue
            sweep_keys.append(_SweepKey(segment, pt))

        sweep_keys = list(sorted(sweep_keys))

        # Add new events for neighbours
        if len(sweep_keys) == 0:
            # If we just removed stuff, the neighbours might now meet...
            if neighbours[0] is not None and neighbours[1] is not None:
                ipt = _nonparallel_intersection_point(
                    neighbours[0].segment, neighbours[1].segment
                )
                if ipt and ipt[1] > pt[1]:
                    add_event(ipt)

            continue

        if neighbours[0] is not None:
            ipt = _nonparallel_intersection_point(
                sweep_keys[0].segment, neighbours[0].segment
            )
            # hyp.note(fstr('IPTL', ipt, pt))
            if ipt and ipt[1] > pt[1]:
                add_event(ipt)

        if neighbours[1] is not None:
            ipt = _nonparallel_intersection_point(
                sweep_keys[-1].segment, neighbours[1].segment
            )
            # hyp.note(fstr('IPTR', ipt, pt))
            if ipt and ipt[1] > pt[1]:
                add_event(ipt)

        # Add them in and continue
        for sweep_key in sweep_keys:
            active.add(sweep_key)


def _det(va, vb):
    return va[0] * vb[1] - va[1] * vb[0]


def _nonparallel_intersection_point(a, b):
    """
    Returns the intersection point or None
    if the segments are parallel or do not intersect
    """
    da = (a[1][0] - a[0][0], a[1][1] - a[0][1])
    db = (b[1][0] - b[0][0], b[1][1] - b[0][1])

    denom = _det(da, db)
    if denom == 0:
        return None

    dab0 = (a[0][0] - b[0][0], a[0][1] - b[0][1])

    t_num = _det(db, dab0)
    u_num = _det(da, dab0)

    if denom < 0:
        t_num = -t_num
        u_num = -u_num
        denom = -denom

    if 0 <= t_num <= denom and 0 <= u_num <= denom:
        ta = t_num / denom
        inter = (a[0][0] + ta * da[0], a[0][1] + ta * da[1])
        return inter

    return None


def _assert_fully_sorted(alist, y):
    """Check that the active list is fully ordered
    """
    for i in range(len(alist)):
        ys = tuple(sorted((alist[i].segment[0][1], alist[i].segment[1][1])))
        assert ys[0] <= y <= ys[1], (ys, y, alist[i])

    for i in range(len(alist)):

        assert alist[i] == alist[i]

        if i > 1:
            assert alist[i - 1] < alist[i]

        # Check that all comparisons are correct
        for j in range(i):
            cmps = (
                alist[i] != alist[j],
                alist[i] == alist[j],
                alist[i] > alist[j],
                alist[i] < alist[j],
                alist[j] != alist[i],
                alist[j] == alist[i],
                alist[j] > alist[i],
                alist[j] < alist[i],
            )
            assert cmps == (
                True,
                False,
                True,
                False,
                True,
                False,
                False,
                True,
            ), (i, j, cmps)
