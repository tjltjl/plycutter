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

"""Representation of the gradually made choices about a plycutter model.
"""

import pyrsistent as pyr
import logging

from .misc import Autocompletable
from .plytypes import fastr, F, Fraction
from .geometry.geom1d import Geom1D

logger = logging.getLogger(__name__)


#
check_interval = 150
check_index = 0


class SheetBuild(pyr.PRecord, Autocompletable):
    """Data evolved by heuristics

    Semantics:

        chosen = set in stone

        ok = as far as we know, this could be

    chosen <= ok always.

    When chosen == ok, we are done

    Attributes

    sheet_chosen -- PMap(sheet_id, Geom2D). Sheet areas chosen
        to be included.
    sheet_ok -- PMap(sheet_id, Geom2D). Sheet areas that are
        not forbidden by e.g.  other sheets or by not being included at all.

    interside_chosen -- PMap(sheet_id, Geom1D). For each interside, the regions
        where the interside's sheet has been chosen.

    interside_ok -- PMap(sheet_id, Geom1D). For each interside, the regions
        where it would still be ok to choose that side.

    There is a variety of sanity checks that can be made on a sheetbuild.
    The checks are relatively slow so by default they are made
    only infrequently (see sheetbuild.check_interval).

    Call sb.check(force=True) to force a check to happen.
    """
    sheetplex = pyr.field()

    sheet_chosen = pyr.field()
    sheet_ok = pyr.field()

    interside_chosen = pyr.field()
    interside_ok = pyr.field()

    dbg = pyr.field()

    def unchoose_interside(self, interside_id, v):
        """State that the range v will never be chosen"""
        assert (self.interside_chosen[interside_id] & v).is_empty(), fastr(
            interside_id, v, (self.interside_chosen[interside_id] & v).area()
        )
        interside = self.sheetplex.interside(interside_id)

        return SheetBuild(
            sheetplex=self.sheetplex,
            sheet_chosen=self.sheet_chosen,
            sheet_ok=self.sheet_ok,
            interside_chosen=self.interside_chosen,
            interside_ok=map_sub_geom(self.interside_ok,
                                      interside.id, v),
        )

    def unchoose_sheet_area(self, sheet_id, area):
        assert (self.sheet_chosen[sheet_id] & area).area() < 1e-5
        return self.transform(['sheet_ok', sheet_id],
                              lambda orig: orig - area)

    def unchoose_intersides(self, unchoices):
        for interside_id, v in unchoices.items():
            self = self.unchoose_interside(interside_id, v)
        return self

    def choose_interside(self, interside_id, v, tentative=False):
        assert (~self.interside_ok[interside_id] & v).measure1d() < 1e-6, \
            (interside_id, self.interside_ok[interside_id], v,
             'INTER',
             (~self.interside_ok[interside_id] & v)
             )
        interside = self.sheetplex.interside(interside_id)
        opposite = interside.opposite(self.sheetplex)

        opposite_2d = opposite.make_fingers(v)
        assert ((opposite_2d &
                 self.sheet_chosen[opposite.sheet]).area() < 1e-6),\
            ((opposite_2d & self.sheet_chosen[opposite.sheet]).area(),
             interside.sheet, opposite.sheet,
             )

        return SheetBuild(
            sheetplex=self.sheetplex,
            sheet_chosen=self.sheet_chosen if tentative else
            # Can't absolutely make the choice
            # due to intersections
            map_or_geom(
                self.sheet_chosen,
                interside.sheet,
                interside.make_fingers(v) & self.sheet_ok[interside.sheet]),
            sheet_ok=map_sub_geom(
                self.sheet_ok,
                opposite.sheet,
                opposite.make_fingers(v)),
            interside_chosen=map_or_geom(
                self.interside_chosen, interside_id, v),
            interside_ok=map_sub_geom(self.interside_ok, opposite.id, v),
        )

    def choose_intersides(self, choices, tentative=False):
        for interside_id, v in choices.items():
            self = self.choose_interside(interside_id, v, tentative=tentative)
        return self

    def check(self, force=False):
        """Assert that various invariants are true.
        """
        global check_index
        check_index += 1
        if force or check_index > check_interval:
            check_index = 0
            self._check()
        else:
            logger.info('check elided')

    def _check(self):
        logger.info('check really')
        sp = self.sheetplex

        # Check that sheet_chosen <= sheet_ok
        for sheet in sp.sheets.values():
            # The slack is allowed but should be removed now that we have
            # exact geometery
            assert (self.sheet_chosen[sheet.id]
                    - self.sheet_ok[sheet.id]).area() < 1e-7, \
                (self.sheet_chosen[sheet.id]
                 - self.sheet_ok[sheet.id]).area()

        for interside in sp.intersides():
            opposite = interside.opposite(sp)
            # Check that interside_chosen <= interside_ok
            assert (self.interside_chosen[interside.id]
                    - self.interside_ok[interside.id]).measure1d() < 1e-7, \
                (interside.id,
                 (self.interside_chosen[interside.id]
                  - self.interside_ok[interside.id]).measure1d(),
                 self.interside_chosen[interside.id],
                 self.interside_ok[interside.id],
                 )

            # Check that intersection can be chosen only by one side
            assert (self.interside_chosen[interside.id] &
                    self.interside_ok[opposite.id]).is_empty()

            # Check that interside_ok for an opposing side is unable
            # to choose parts that have been chosen in the sheet.
            chosen_but_other_ok = (self.sheet_chosen[interside.sheet] &
                                   interside.make_fingers(
                self.interside_ok[opposite.id]))
            assert chosen_but_other_ok.area() < 0.001, \
                fastr(interside.id,
                      self.interside_ok[interside.id],
                      self.interside_ok[opposite.id],
                      self.interside_chosen[interside.id],
                      self.interside_chosen[opposite.id],
                      chosen_but_other_ok.area())

        logger.info('check done')
        return self


def map_or_geom(mapping, key, geom):
    """For a pyrsistent map, result of map[key] |= geom"""
    return mapping.set(key,
                       (mapping.get(key) or geom.__class__.empty()) | geom)


def map_and_geom(mapping, key, geom):
    """For a pyrsistent map, result of map[key] &= geom"""
    return mapping.set(key, mapping.get(key) & geom)


def map_sub_geom(mapping, key, geom):
    """For a pyrsistent map, result of map[key] -= geom"""
    return mapping.set(key, mapping.get(key) - geom)


def map_and(mapping1, mapping2):
    """Return the intersection of two arbitrary maps whose values are geoms.
    """
    return pyr.pmap({
        k: v1 & mapping2[k]
        for k, v1 in mapping1.items()
        if k in mapping2 and not (v1 & mapping2[k]).is_empty()
    })


def map_or(mapping1, mapping2):
    """Return the union of two arbitrary maps whose values are geoms.
    """
    for k, v in mapping2.items():
        mapping1 = map_or_geom(mapping1, k, v)
    return mapping1


def create_sheetbuild(sp, params):
    """Create an as-free-as-possible SheetBuild for the given SheetPlex.
    """
    max_buffer = 0  # F(1, 1000)
    both_faces_buffer = 1
    sheet_chosen = {sheet.id: (sheet.faces[0] & sheet.faces[1]).buffer(
        both_faces_buffer) & sheet.slices_max.buffer(max_buffer)
        for sheet in sp.sheets.values()}
    sheet_chosen_orig = {**sheet_chosen}

    ok_rad = params['support_radius']
    sheet_ok = {sheet.id:
                sheet.faces[0].buffer(ok_rad)
                & sheet.faces[1].buffer(ok_rad)
                & sheet.slices_max.buffer(max_buffer)
                for sheet in sp.sheets.values()}

    interside_chosen = {interside.id: Geom1D.empty()
                        for interside in sp.intersides()}

    interside_ok = {interside.id:
                    interside.project_to_1d(
                        sheet_ok[interside.sheet] &
                        sp.sheets[interside.sheet].slices_max &
                        interside.make_fingers(
                            Geom1D([[Fraction(-1000), Fraction(1000)]])))
                    for interside in sp.intersides()}

    for interside in sp.intersides():
        interside_ok[interside.id] &= interside_ok[interside.opposite(sp).id]

    for k, v in interside_ok.items():
        interside = sp.interside(k)
        sheet_chosen[interside.sheet] -= interside.make_fingers(v)

    global locs
    locs = locals()

    return SheetBuild(
        sheetplex=sp,
        sheet_chosen=pyr.pmap(sheet_chosen),
        sheet_ok=pyr.pmap(sheet_ok),
        interside_chosen=pyr.pmap(interside_chosen),
        interside_ok=pyr.pmap(interside_ok),
    )


def show_sheet(ax, sheetplex, sheetbuild, sheet_id):
    """Show a sheet build progress in a matplotlib axes object.
    """
    sheetbuild.sheet_ok[sheet_id].show2d(ax, 'green', alpha=1.0)
    sheetbuild.sheet_chosen[sheet_id].show2d(ax, 'cyan', linewidth=1)

    for interside in sheetplex.intersides(sheet_id):
        interside.make_fingers(sheetbuild.interside_chosen[interside.id],
                               F(4, 10), F(6, 10)).\
            show2d(ax, 'red', alpha=0.3)
        interside.make_fingers(sheetbuild.interside_ok[interside.id],
                               F(3, 10), F(7, 10)).\
            show2d(ax, 'blue', alpha=0.3)
