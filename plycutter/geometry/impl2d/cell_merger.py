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

import collections
from itertools import islice

from ...plytypes import fstr


def ifirst(iterator):
    return list(islice(iterator, 1))[0]


def canonicalize_cell(cell):
    mini = 0
    for i in range(len(cell)):
        if cell[i] < cell[mini]:
            mini = i
    return list(cell[mini:]) + list(cell[:mini])


class CellMergerSlow:
    """Given directed circular vertex lists, merge opposite half-edges.

    Keeps track of the original cells merged into a given
    final cell.

    """
    class Cell:
        def __init__(self, cell, origs):
            self.vertices = cell
            self.half_edges = collections.defaultdict(lambda: 0)
            self.origs = origs

            n = len(cell)
            for i in range(n):
                v0 = cell[i]
                v1 = cell[(i + 1) % n]
                self.half_edges[v0, v1] += 1

        def __repr__(self):
            return f'CELL({fstr(self.vertices), fstr(self.origs)})'

        def merged(self, other):
            self_half_edges = {**self.half_edges}
            other_half_edges = {**other.half_edges}
            changed = False
            for (v0, v1), v in list(self_half_edges.items()):
                assert (v0, v1) not in other_half_edges
                if (v1, v0) in other_half_edges:
                    del self_half_edges[v0, v1]
                    del other_half_edges[v1, v0]
                    changed = True
            if not changed:
                return None

            # Can have multiple contacts...
            nxt = collections.defaultdict(lambda: set())
            for hes in (self_half_edges, other_half_edges):
                for v0, v1 in hes:
                    nxt[v0].add(v1)

            res = []
            # print('merge', fstr(self_half_edges, other_half_edges))
            # print('nxt', fstr(nxt))
            while len(nxt):
                lst = []
                k0, vs = nxt.popitem()
                lst.append(k0)
                v = ifirst(vs)
                if len(vs) > 1:
                    nxt[k0] = vs - {v}
                while v != k0:
                    # print(v)
                    lst.append(v)
                    vs = nxt.pop(v)
                    vnext = ifirst(vs)
                    if len(vs) > 1:
                        nxt[v] = vs - {vnext}
                    v = vnext
                res.append(CellMergerSlow.Cell(tuple(lst),
                                               self.origs + other.origs))
            return res

    def __init__(self):
        self.cells = set()

        # All the original cells added, for debugging
        self.original_cells = set()
        self.original_cells_list = []

        self.annihilated_original_cells = set()

    def add_cell(self, cell):
        assert cell not in self.original_cells
        self.original_cells.add(cell)
        self.original_cells_list.append(cell)

        cell = self.Cell(cell, (cell,))
        # print('Add', fstr(cell))

        to_add = [cell]
        while len(to_add):
            cell = to_add.pop()

            deleted = False

            for other in list(self.cells):
                m = cell.merged(other)
                if m is not None:
                    self.cells.remove(other)
                    if len(m) == 0:
                        # Special case: annihilation.
                        # Need to add both original cells
                        # to all other groups with those original
                        # cells

                        assert len(set(cell.origs) & set(other.origs)) == 0
                        orig_cells = set(cell.origs) | set(other.origs)
                        saved = False
                        for outside_cell in self.cells:
                            if set(outside_cell.origs) & orig_cells:
                                # Changes the value inside the obj
                                # --> the list above is ok
                                outside_cell.origs = tuple(
                                    set(outside_cell.origs) | orig_cells)
                                saved = True

                        if not saved:
                            self.annihilated_original_cells |= orig_cells
                        deleted = True
                        break
                    cell = m[0]
                    to_add.extend(m[1:])

            if deleted:
                continue

            self.cells.add(cell)

        # print('Cells after', fstr(self.cells))

    def get_cells(self):
        self.check()
        return [cell.vertices for cell in self.cells]

    def get_cells_originals(self):
        self.check()
        return {cell.vertices: cell.origs for cell in self.cells}

    def check(self):
        origs = set()
        for cell in self.cells:
            origs |= set(cell.origs)
        assert not (origs & self.annihilated_original_cells)
        origs |= self.annihilated_original_cells

        if origs != self.original_cells:

            assert origs == self.original_cells, [
                self.original_cells_list,  (self.original_cells, origs)]
