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
import hypothesis.stateful
import hypothesis.strategies as hys

from ..segment_tree import SegmentTree1D, SegmentTree1D_Slow


@hyp.given(
    hys.lists(elements=hys.tuples(
        # hys.fractions(),
        # hys.fractions(),
        hys.integers(),
        hys.integers(),
    )),
    hys.lists(elements=hys.fractions()))
def test_segment_tree_1d(segments, points):
    i = [-1]

    def mkid():
        i[0] += 1
        return f"s{i[0]}"

    segments = [
        tuple(sorted((s[0], s[1]))) + (mkid(),)
        for s in segments
    ]

    st_slow = SegmentTree1D_Slow(segments)
    st = SegmentTree1D(segments)

    def pts():
        yield from iter(points)
        for low, high, id in segments:
            yield low
            yield high

    for pt in pts():
        v0 = list(st_slow.stab(pt))
        v1 = list(st.stab(pt))

        assert list(sorted(v0)) == list(sorted(v1)), (pt, str(st.root))
