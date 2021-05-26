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

import logging
import numpy as np
import operator
import math
from itertools import islice
from sortedcontainers import SortedDict
import functools

from ...plytypes import Fraction, FractionClass, FractionFromExact, \
    FractionFromFloat, fstr, totuple
from ..geom1d import Geom1D

from .primitive2d import edges_iter, simple_polygon_area, \
    locate_point_polygon_winding_number
from .. import geom2dbase
from . import arrangement2d
from .frangle2d import frangle_unit_square_vector, vector_frangle, \
        MAX_FRANGLE, FRANGLE_180


def NP(v): return np.array(v, object)


def L(v): return np.array(v, object).tolist()


def buffer_of_size(sz, n=8):
    sz = FractionFromExact(sz)
    buffer_polygon = []
    for i in range(n):
        # Use frangles to be exactly
        # the same and symmetric on all quadrants
        frangle = i * Fraction(MAX_FRANGLE) / n
        v = frangle_unit_square_vector(frangle)
        vlen2 = sum(NP(v) ** 2)
        # This is the approximate vector length.
        # TODO: search for pythagorean triples to
        # avoid this :)
        vlen = Fraction(FractionFromFloat(
                math.sqrt(float(vlen2)), -10)) / sz
        buffer_polygon.append((v[0] / vlen, v[1] / vlen))
    return buffer_polygon


def fr_points(points):
    res = []
    for p in points:
        assert len(p) == 2
        res.append((Fraction(p[0]), Fraction(p[1])))
    return res


def check_no_dups(verts):
    for v0, v1 in edges_iter(verts):
        assert np.any(v0 != v1), (v0, v1, verts)


class Simple_Polygon_With_Holes:
    """Requires clockwise polygons
    """

    def __init__(self, outer, holes=[]):
        check_no_dups(outer)
        a = simple_polygon_area(outer)
        assert a > 0, a
        self.outer = outer
        for hole in holes:
            assert simple_polygon_area(hole) > 0
            check_no_dups(hole)
        self.holes = holes

    def repro(self):
        """Gets so long that not default"""
        return "plycutter.geom2d_py.Simple_Polygon_With_Holes(%s,%s)" % (
            L(self.outer), L(self.holes)
        )

    def transformed_with(self, f):
        def fpoly(p):
            p = f(p)
            if simple_polygon_area(p) < 0:
                p = list(reversed(p))
            return p
        return self.__class__(
            fpoly(self.outer),
            [fpoly(hole) for hole in self.holes])

    def segments(self):
        yield from edges_iter(self.outer)
        for hole in self.holes:
            yield from edges_iter(hole)

    def area(self):
        a = simple_polygon_area(self.outer)
        for hole in self.holes:
            a = a - simple_polygon_area(hole)
        return a

    def locate(self, point):
        w = locate_point_polygon_winding_number(self.outer, point)
        # print('loc', fstr(point, w, approx=True))
        if w == -1:
            return w
        if w == 0:
            return w

        for hole in self.holes:
            h = locate_point_polygon_winding_number(hole, point)
            # print('hloc', fstr(point, h, approx=True))
            if h == 0:
                return 0
            if h > 0:
                return -1
        return 1


oplogger = logging.getLogger("geom2dops")
oplogger.setLevel(logging.INFO)

arrangement_check = True

UNION_N = 64


class Geom2D_Py_Simple_Holes(geom2dbase.Geom2DBase):
    def __init__(self, spwhs):
        assert len(spwhs) >= 0  # "typecheck"
        self.spwhs = spwhs

    def repro(self):
        """This gets so long that it is not the default"""
        return "plycutter.geom2d_py.Geom2D_Py_Simple_Holes([%s])" % \
            ",".join([spwh.repro() for spwh in self.spwhs])

    def _arrangement(self, geoms):
        global ggeoms
        ggeoms = geoms
        segs = {}
        n = 0
        for geom in geoms:
            for segment in geom.all_segments():
                segs['s%d' % n] = segment
                n += 1

        arr = arrangement2d.SegmentArrangement(segs)
        return arr

    def _arrangement_op(self, geoms, op):
        arr = self._arrangement(geoms)

        cells = set()
        cell_points = arr.cell_points()
        for cell, pt in cell_points.items():
            locs = [geom.locate(pt) for geom in geoms]
            for loc in locs:
                assert loc == -1 or loc == 1
            # print('ARRCELL', fstr(cell), locs)
            if op(*[loc > 0 for loc in locs]):
                cells.add(cell)

        polys = arr.poly_cells(cells)
        result = Geom2D_Py_Simple_Holes(polys)

        if arrangement_check:
            # Check that we really got what we nitended
            for cell, pt in cell_points.items():
                locs = [geom.locate(pt) for geom in geoms]
                for loc in locs:
                    assert loc == -1 or loc == 1
                # print('ARRCELL', fstr(cell), locs)
                new_loc = result.locate(pt)
                if op(*[loc > 0 for loc in locs]):
                    assert new_loc > 0, fstr(pt, new_loc, approx=True)
                else:
                    assert new_loc < 0, fstr(pt, new_loc, approx=True)

        return result

    @classmethod
    def empty(cls):
        return Geom2D_Py_Simple_Holes(
            [])

    @classmethod
    def union(cls, geoms):
        if len(geoms) > UNION_N:
            return union_tree(geoms)
        return cls.empty()._arrangement_op(geoms,
                                           lambda *args: any(args))

    def is_empty(self):
        return len(self.spwhs) == 0

    def area(self):
        return sum(spwh.area() for spwh in self.spwhs)

    @classmethod
    def polygon(cls, points, reorient=False):
        # print(points)
        for point in points:
            for coord in point:
                assert isinstance(coord, FractionClass) or int(
                    coord) == coord, (coord, points)
        points = fr_points(points)
        # print(points)

        if reorient and simple_polygon_area(points) < 0:
            points = list(reversed(points))

        result = Geom2D_Py_Simple_Holes(
            [Simple_Polygon_With_Holes(points)])
        assert result.area() >= 0, result.area()
        return result

    @classmethod
    def from_shapely(cls, sh):
        if sh.geom_type == 'GeometryCollection':
            geoms = [cls.from_shapely_polygon(geom)
                     for geom in sh.geoms]
            spwhs = []
            for geom in geoms:
                assert len(geom.spwhs) == 1
                spwhs.extend(geom.spwhs)
            return cls(spwhs)
        if sh.geom_type == 'Polygon':
            return cls.from_shapely_polygon(sh)
        if sh.geom_type == 'MultiPolygon':
            geoms = [cls.from_shapely_polygon(poly)
                     for poly in sh.geoms]
            spwhs = []
            for geom in geoms:
                assert len(geom.spwhs) == 1
                spwhs.extend(geom.spwhs)
            return cls(spwhs)
        raise Exception(f'Invalid geom type {sh.geom_type}')

    @classmethod
    def from_shapely_polygon(cls, sh):
        import shapely
        sh = shapely.geometry.polygon.orient(sh, -1)

        def sh_points(pts):
            pts = fr_points(pts)
            assert pts[0] == pts[-1]
            return pts[:-1]
        result = cls([Simple_Polygon_With_Holes(
            sh_points(sh.exterior.coords),
            [list(reversed(sh_points(hole.coords)))
             for hole in sh.interiors])])
        assert result.area() >= 0, result.area()
        return result

    def __or__(self, other):
        if self.is_empty():
            if oplogger.isEnabledFor(logging.DEBUG):
                oplogger.debug('or shortcut %s', other.area())
            return other
        if other.is_empty():
            if oplogger.isEnabledFor(logging.DEBUG):
                oplogger.debug('or shortcut %s', self.area())
            return self
        if oplogger.isEnabledFor(logging.DEBUG):
            oplogger.debug('or %s %s', self.area(), other.area())
        return self._arrangement_op([self, other],
                                    operator.or_)

    def __and__(self, other):
        if oplogger.isEnabledFor(logging.DEBUG):
            oplogger.debug('and %s %s', self.area(), other.area())
        return self._arrangement_op([self, other],
                                    operator.and_)

    def __sub__(self, other):
        if other.is_empty():
            if oplogger.isEnabledFor(logging.DEBUG):
                oplogger.debug('sub shortcut %s', self.area())
            return self
        if oplogger.isEnabledFor(logging.DEBUG):
            oplogger.debug('sub %s %s', self.area(), other.area())

        def oper(a, b):
            # print('suboper', a, b)
            if a and not b:
                return True
            return False

        return self._arrangement_op([self, other], oper)

    def buffer(self, amount, resolution=8):
        global locs
        locs = locals()  # DBG
        if oplogger.isEnabledFor(logging.DEBUG):
            oplogger.debug('buffer %s %s', amount,
                           len(list(self.all_segments())))
        if amount == 0:
            return self
        amount = FractionFromExact(amount)
        buf = buffer_of_size(abs(amount), n=resolution)
        res = self._minkowski_convex_op(
            buf, 1 if amount > 0 else -1)
        oplogger.debug('done')
        return res

    def all_segments(self):
        for poly in self.spwhs:
            yield from poly.segments()

    def locate(self, pt):
        if len(self.spwhs) == 0:
            return -1

        return max([poly.locate(pt) for poly in self.spwhs])

    def polygons(self):
        return [self.__class__([spwh]) for spwh in self.spwhs]

    def exterior(self):
        assert len(self.spwhs) == 1
        return self.spwhs[0].outer

    def holes(self):
        assert len(self.spwhs) == 1
        yield from self.spwhs[0].holes

    def project_to_1d(self, vec, offs):
        res = Geom1D.empty()
        for spwh in self.spwhs:
            ds = np.dot(np.array(spwh.outer, object), vec) + offs
            res = res | Geom1D([(min(ds), max(ds))])
        return res

    def transformed_with(self, f):
        return self.__class__(
            [spwh.transformed_with(f) for spwh in self.spwhs])

    def _minkowski_convex_op(self, cpoly, sign):
        cpoly = NP(cpoly)
        assert sign != 0
#        if sign < 0:
#            cpoly = cpoly * -1
        offset = NP((0, 0))
        if locate_point_polygon_winding_number(cpoly, offset) != 1:
            offset = sign * NP(cpoly[0] + cpoly[1] + cpoly[2]) / 3
            cpoly = cpoly - sign * offset
            assert locate_point_polygon_winding_number(cpoly, (0, 0)) == 1
        # print('offset', offset)

        pos = []
        neg = []

        holes = []

        for spwh in self.spwhs:
            pos.append(Geom2D_Py_Simple_Holes(
                [Simple_Polygon_With_Holes(
                    spwh.transformed_with(lambda v: v + offset).outer)]))

            ext_simples = minkowski_ribbon_simple_convex(
                spwh.outer, cpoly, sign)
            holes.extend(spwh.holes)

            def reorient(e):
                if simple_polygon_area(e) < 0:
                    return list(reversed(e))
                return e
            ext = [Geom2D_Py_Simple_Holes(
                [Simple_Polygon_With_Holes(reorient(e) + offset)])
                for e in ext_simples]
            if sign > 0:
                pos.extend(ext)
            else:
                neg.extend(ext)

        for hole in holes:
            g = Geom2D_Py_Simple_Holes([Simple_Polygon_With_Holes(hole)])
            g = g._minkowski_convex_op(cpoly, -sign)
            neg.append(g)
            # print('Hole', g.area())

        if True:
            # Reveals weird bugs...
            gpos = union_tree(pos)
            gneg = union_tree(neg)
            geom = gpos - gneg

            return geom

        npos = len(pos)
        nneg = len(neg)

        posmask = np.arange(npos + nneg) < npos
        negmask = ~posmask

        subs = pos + neg

        return self._arrangement_op(
                subs,
                lambda *v: np.any(posmask & v) & ~np.any(negmask & v))


Geom2D_Py_Simple_Holes.__doc__ = geom2dbase.Geom2DBase.__doc__

error_geoms = None


def union_tree(geoms):
    # N = 8
    while len(geoms) > 1:
        # print(sum([g.area() for g in geoms]), len(geoms))
        nxt = []
        for i in range(0, len(geoms), UNION_N):
            slc = geoms[i:i + UNION_N]
            if len(slc) > 1:
                try:
                    nxt.append(Geom2D_Py_Simple_Holes.union(slc))
                except Exception:
                    global error_geoms
                    error_geoms = slc
                    raise
            else:
                assert len(slc) == 1
                nxt.append(slc[0])
        # print('A', len(nxt))
        geoms = nxt
    # print(sum([g.area() for g in geoms]), len(geoms))
    if len(geoms) == 0:
        return Geom2D_Py_Simple_Holes.empty()
    return geoms[0]


def _ifirst(iterator):
    return list(islice(iterator, 1))[0]


def _triple_vertices_iter(cell):
    n = len(cell)
    for i in range(n):
        yield (cell[(i - 1) % n], cell[i], cell[(i + 1) % n])

def minkowski_ribbon_simple_convex(spoly, cpoly, sign):  # noqa # XXX
    """Generate the simple polygons, the union of which
    are added or semoved to spoly.

    Assuming clockwise polygon that contains the origin
    """
    # XXX This is way too complex...
    assert len(cpoly) > 2
    assert sign == -1 or sign == 1

    # print(fstr(spoly))
    # print(fstr(cpoly))
    # print(sign)

    edge_start_by_frangle = SortedDict()
    for i, (v0, v1) in enumerate(edges_iter(cpoly)):
        frangle = vector_frangle(NP(v1) - v0)
        assert frangle >= 0
        if frangle in edge_start_by_frangle:
            # Leave out if it is already in (degenerate
            # vertex at 180 degree angle)
            continue
        edge_start_by_frangle[frangle] = i
        # Avoid special cases by rolling this over
        edge_start_by_frangle[frangle + MAX_FRANGLE] = i
        edge_start_by_frangle[frangle - MAX_FRANGLE] = i

    parts = []

    def find_vertex_on(frangle):
        fr = _ifirst(edge_start_by_frangle.irange(
            -8, frangle, inclusive=(True, True), reverse=True))
        return (edge_start_by_frangle[fr] + 1) % len(cpoly)

    # Generate the corner pieces + v..vn edge piece
    for vp, v, vn in _triple_vertices_iter(spoly):
        pfr = vector_frangle(NP(v) - vp)
        nfr = vector_frangle(NP(vn) - v)

        # Find the cpoly vertex visible on the edges
        iprev = find_vertex_on(pfr)
        inext = find_vertex_on(nfr)

        others = [NP(vn) + sign * NP(cpoly[inext]), vn, v]
        if (nfr - pfr) % MAX_FRANGLE < FRANGLE_180:
            # Turn right
            idxs = np.arange(iprev, 1 + (inext if inext >= iprev
                                         else inext + len(cpoly)))
            if sign < 0:
                idxs = np.array([inext])
        else:
            # Turn left
            idxs = np.arange(inext, 1 + (iprev if iprev >= inext
                                         else iprev + len(cpoly)))
            if sign > 0:
                idxs = np.array([inext])
            others = list(reversed(others))

        idxs = idxs % len(cpoly)

        # print(fstr('VERT', (vp, v, vn)))
        # print(fstr(iprev, inext, idxs, others))

        if sign < 0:
            idxs = idxs[::-1]
            others = list(reversed(others))

        if len(idxs) == 0:
            continue

        # print(idxs)

        verts = [NP(v) + sign * NP(cpoly[idx]) for idx in idxs]
        # If non-convex by being too sharp corner (i.e. vertex would
        # "poke through", add the vertex there
        adds = set()
        for i in range(0, len(verts) - 1):
            if sign * simple_polygon_area((verts[i], verts[i + 1], v)) < 0:
                adds.add(i)
        if len(adds) > 0:
            # print('Reverse', sorted(adds))
            for i in reversed(sorted(adds)):
                verts[i:i] = [v]

        verts = verts + others
        # print(fstr(verts))
        for vert in verts:
            assert len(vert) == 2, vert

        verts = simplify_degeneracies(verts)
        # print(fstr('SIMP', verts))

        if len(verts) > 2:
            # Generate the full polygon
            parts.append(verts)
        else:
            pass
            # print('NOGEN')

    return parts


def simplify_degeneracies(verts):
    changed = True
    verts = list(totuple(verts))

    @functools.lru_cache(maxsize=None)
    def area(a, b, c):
        return simple_polygon_area([a, b, c])

    while changed:
        changed = False

        # Whenever the area of three vertices is zero,
        # we can remove the middle one

        i = 0
        nz = False
        while i < len(verts):
            n = len(verts)
            vpp = verts[(i - 2) % n]
            vp = verts[(i - 1) % n]
            v = verts[i]
            vn = verts[(i + 1) % n]

            if area(vp, v, vn) != 0:
                nz = True

            # print(fstr(i, v, area(vpp, vp, v), area(vp, v, vn)))
            if area(vpp, vp, v) != 0 and area(vp, v, vn) == 0:
                # Delete
                verts[i:i+1] = []
                changed = True
            else:
                i += 1

        # If all areas are zero, remove first.
        # Could still be a real polygon with just all vertices repeated
        if not nz and len(verts) > 0:
            verts[0:1] = []
            changed = True

    return verts
