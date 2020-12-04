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

from ..cell_merger import CellMergerSlow, canonicalize_cell
import numpy as np
import hypothesis as hyp
import hypothesis.strategies as hys

from ....plytypes import Fraction, totuple

F = Fraction


hyp.settings.register_profile('verbose', verbosity=hyp.Verbosity.verbose,
                              max_examples=1000)
hyp.settings.load_profile('verbose')

SQUARE = np.array([
    [0, 0],
    [0, 1],
    [1, 1],
    [1, 0],
])

N = 3
SQUARES = {
    (x, y): SQUARE + [x, y]
    for x in range(N)
    for y in range(N)
}

SQUARES_KEYS = list(sorted(SQUARES.keys()))


@hyp.settings(deadline=30000)
@hyp.given(hys.lists(min_size=2, max_size=20,
                     elements=hys.permutations(SQUARES_KEYS)))
#        elements=hys.sampled_from(SQUARES_KEYS), unique=True))))
def test_cell_merger(perms):
    canonical = None
    for perm in perms:
        m = CellMergerSlow()
        for k in perm:
            cell = totuple(SQUARES[k])
            m.add_cell(cell)
            hyp.note(('add', cell))

            if False:
                origs_used = set()
                for k, v in m.get_cells_originals().items():
                    hyp.note(k)
                    for orig in v:
                        hyp.note('    %s' % (orig,))
                    vset = set(v)
                    assert len(origs_used & vset) == 0
                    origs_used |= vset

        cells = list(sorted([canonicalize_cell(c) for c in m.get_cells()]))

        if canonical is None:
            canonical = cells

        assert cells == canonical

        assert len(cells) == 1
