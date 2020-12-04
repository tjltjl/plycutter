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

"""Representation for sheets extracted from a 3D model.
"""

import pyrsistent as pyr
import numpy as np
import logging

from .plytypes import FractionClass, FractionFromExact
from .misc import Autocompletable
from .geometry.geom2d import Geom2D

logger = logging.getLogger(__name__)


class SheetPlex(pyr.PRecord, Autocompletable):
    """Constant data extracted from a 3D model.

    Three classes represent the objects that live in a sheetplex:

    * `Sheet` -- a single plywood (or other material) sheet
    * `Intersection` -- the intersection between two sheets
    * `InterSide` -- one sheet's side of an Intersection

    The direct members are

        intersections -- PMap[str, Intersection]

    Stored in a pyrsistent.PRecord, see the methods
    therein for how to manipulate this object.
    """

    sheets = pyr.field(type=pyr.PMap)
    """A PMap of all sheets by their id."""

    intersections = pyr.field(type=pyr.PMap)
    """A PMap of all intersections by their id."""

    def interside(self, id):
        """Get the interside with the given id"""
        inter_id, sheet_id = id
        sides = self.intersections[inter_id].sides
        if sides[0].sheet == sheet_id:
            return sides[0]
        else:
            assert sides[1].sheet == sheet_id
            return sides[1]

    def intersides(self, sheet_id=None):
        """Return an iterator for intersides of a particular sheet or all"""
        for intersection in self.intersections.values():
            for side in intersection.sides:
                if sheet_id is None or side.sheet == sheet_id:
                    yield side

    def limit(self, sheet_ids):
        """Return a sub-Sheetplex with only the given sheet_ids included"""
        sheet_ids = set(sheet_ids)
        return SheetPlex(
            sheets=pyr.pmap({
                sheet_id: sheet for sheet_id, sheet in self.sheets.items()
                if sheet_id in sheet_ids
            }),
            intersections=pyr.pmap({
                inter_id: inter
                for inter_id, inter
                in self.intersections.items()
                if inter.sides[0].sheet in sheet_ids
                and iter.sides[1].sheet in sheet_ids
            })
        )


class Sheet(pyr.PRecord, Autocompletable):
    """A single sheet that is part of the full model.

    Contains various data extracted from the 3D model
    to help geometric reasoning about the joints.
    """

    # Initial

    id = pyr.field(type=str)
    sides = pyr.field()
    slice_rotation = pyr.field()
    inverse_slice_rotation = pyr.field()

    # 2D

    slices = pyr.field()
    slices_height = pyr.field()
    faces = pyr.field()
    outside_slices = pyr.field()

    # Derived

    slices_max = pyr.field()
    slices_max = pyr.field()
    slices_min = pyr.field()
    both_faces = pyr.field()

    def normal(self, side):
        return self.sides[side][0:3]

    def offset(self, side):
        return self.sides[side, 3]

    def midplane(self):
        return 0.5 * np.array([1, -1]).dot(self.sides)

    def project_point4(self, pt):
        """Project a homogeneous point or line to this plane,
        return 2D point or direction vector for a line"""
        p4 = self.inverse_slice_rotation.dot(pt)
        p2 = p4[0:2]
        if p4[3] != 0:
            p2 /= p4[3]
        return p2


class Intersection(pyr.PRecord, Autocompletable):
    """An intersection between two sheets.

    An intersection defines the *intersection coordinate system*
    which is a 1D coordinate across the intersection.

    The joint fingers are defined in terms of that coordinate.
    """
    id = pyr.field(type=str)
    """The identifier of this intersection
    """

    direction = pyr.field()
    """The 3-vector direction of the intersection in 3-space.
    """

    origins = pyr.field()
    """Origins of all intersection lines between the two intersecting sheets.
    """

    sides = pyr.field(type=tuple)
    """The ``InterSide`` objects of this intersection.

    The ``InterSide`` objects store information about the location
    of the intersection within the sheets themselves.
    """

    def sheets(self):
        """Return a tuple of the particpating sheet ids """
        return [interside.sheet for interside in self.sides]

    def side_by_sheet(self, sheet_id):
        """Get the ``InterSide`` for the given sheet_id.

        Returns:
            None if the given ``sheet_id`` is not in this intersection.
        """
        for side in self.sides:
            if side.sheet == sheet_id:
                return side
        return None


def robust_isfinite(v):
    """A version of isfinite that also works on huge rationals"""
    return v + 1 > v


class InterSide(pyr.PRecord, Autocompletable):
    """A particular sheet's view of an intersection.
    """
    id = pyr.field(type=tuple)
    """The identifier of the InterSide.
    """

    intersection = pyr.field(type=str)
    """The id of the ``Intersection`` this ``InterSide`` belongs to
    """

    idx = pyr.field(type=int)
    """The index of this ``InterSide`` within the intersection.
    """

    sheet = pyr.field(type=str)
    """The id of the `Sheet` on which this ``InterSide`` lies.
    """

    direction = pyr.field()
    """The Fraction 2-vector direction on the 2D sheet.
    """

    origin_offset = pyr.field(type=FractionClass)
    """The offset along ``direction`` to the intersections coordinate's origin
    """

    normal = pyr.field()
    """The normal vector to the intersection in the 2D plane.
    """

    min_normal_offset = pyr.field(type=FractionClass)
    """The minimum value of normal dot pt that belongs to the intersection.
    """

    max_normal_offset = pyr.field(type=FractionClass)
    """The maximum value of normal dot pt that belongs to the intersection.
    """

    def opposite(self, sheetplex):
        """Get the other InterSide in the Intersection we are part of"""
        return sheetplex.intersections[self.intersection].sides[1 - self.idx]

    def project_to_1d(self, geom2d):
        """Project a Geom2D in Sheet coords to a Geom1D in Intersection coords.
        """
        return geom2d.project_to_1d(self.direction, -self.origin_offset)

    def make_fingers(self, where, fract0=0, fract1=1):
        """
        Create a geom for fingers on the sheet on this InterSide.

        The parameters fract0 and fract1 are
        useful mostly for visualizations of different rules affecting
        the fingers.

        Parameters:
            where: Geom1D of sections along joint
            fract0: From this fraction coordinate across joint
            fract1: To this fraction coordinate across joint
        Returns:
            2D geometry of fingers.
        """
        res = []
        fract0 = FractionFromExact(fract0)
        fract1 = FractionFromExact(fract1)
        for interval in where.get_intervals():
            assert isinstance(interval[0], FractionClass)
            assert isinstance(interval[1], FractionClass)
            dire = self.direction
            norm = self.normal
            offs = self.origin_offset

            def coords(along, across):
                assert robust_isfinite(along), along
                al = (along + offs)

                fract = fract0 + (fract1 - fract0) * across
                nc = self.min_normal_offset + \
                    (self.max_normal_offset - self.min_normal_offset) * fract
                result = al * dire + nc * norm
                assert np.all([isinstance(r, FractionClass)
                               for r in result]), result
                return result

            # logger.debug('mkfingers %s',
            #    (self.id, fstr(interval, dire, norm, offs, fract0, fract1)))

            res.append(Geom2D.polygon([
                coords(interval[0], 0),
                coords(interval[0], 1),
                coords(interval[1], 1),
                coords(interval[1], 0),
            ], reorient=True))
        return Geom2D.union(res)
