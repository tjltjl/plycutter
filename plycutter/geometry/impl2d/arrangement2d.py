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

import numpy as np
from sortedcontainers import SortedList
from itertools import islice
from typing import Any, Set, Dict
import collections

from ...plytypes import Fraction, FractionClass, fstr, fastr
from ..types2d import Point, Segment
from ..aabb import AABB

from .frangle2d import (
    vector_frangle,
    MAX_FRANGLE,
    FRANGLE_90,
    frangle_unit_square_vector,
)
from .primitive2d import (
    segment_segment_general_intersection,
    line_segment_point_fraction,
    simple_polygon_area,
    locate_point_polygon_winding_number,
    edges_iter,
)
from . import scan2d, segment_tree
from . import cell_merger, pmr_quadtree
from . import geom2d_simple_holes


def NP(v):
    return np.array(v, object)


class _Vertex:
    coordinates: Point
    half_edges_by_order: SortedList

    def _key(self, half_edge):
        assert half_edge[0] == self.coordinates
        return vector_frangle(
            np.array(half_edge[1], dtype=object) - half_edge[0]
        )

    def __init__(self, coordinates):
        self.coordinates = coordinates

        self.half_edges_by_order = SortedList(key=self._key)

    def add_half_edge(self, half_edge):
        existing = list(
            self.half_edges_by_order.irange(
                half_edge, half_edge, inclusive=(True, True)
            )
        )
        assert len(existing) == 0, fstr(half_edge, existing, approx=True)

        self.half_edges_by_order.add(half_edge)

    def get_next_cw(self, half_edge):
        assert half_edge[0] == self.coordinates

        nxt = list(
            islice(
                self.half_edges_by_order.irange(
                    half_edge, None, inclusive=(False, True)
                ),
                1,
            )
        )
        if len(nxt) == 0:
            nxt = list(
                islice(
                    self.half_edges_by_order.irange(
                        None, None, inclusive=(True, True)
                    ),
                    1,
                )
            )

        assert len(nxt) == 1

        return nxt[0]

    def get_next_ccw(self, half_edge):
        assert half_edge[0] == self.coordinates

        nxt = list(
            islice(
                self.half_edges_by_order.irange(
                    None, half_edge, inclusive=(True, False), reverse=True
                ),
                1,
            )
        )
        if len(nxt) == 0:
            nxt = list(
                islice(
                    self.half_edges_by_order.irange(
                        None, None, inclusive=(True, True), reverse=True
                    ),
                    1,
                )
            )

        assert len(nxt) == 1, (self.half_edges_by_order, nxt)

        prev = self.get_next_cw(nxt[0])
        assert prev == half_edge, fstr(
            half_edge,
            nxt[0],
            prev,
            self._key(half_edge),
            self._key(nxt[0]),
            self._key(prev),
        )

        return nxt[0]


class _HalfEdge:
    vs: Segment
    sids: Set
    cell: Any  # Tuple of Points (cell)
    frangle: Fraction

    def __init__(self, vs):
        self.vs = vs
        self.sids = set()
        self.cell = None
        self.frangle = vector_frangle(
            np.array(self.vs[1], object) - self.vs[0]
        )


class _Cell:
    negative: bool
    count_from_outside: int  # Number of cells outside this
    container: Any  # Outside this one
    contains: Any  # For positive cells: what is inside this one
    vertices: Any  # List of vertices

    def __repr__(self):
        return (
            f"_Cell({self.negative}, {self.count_from_outside}, "
            f"container={self.container}, contains={self.contains}, "
            f"verts={self.vertices})"
        )


class SegmentArrangement:
    """Represent a 2D line segment arrangement of full polygons.

    Several opportunities for cleanup and optimization remain.

    Attributes:
        vertices -- all points: dict: id -> _Vertex
        edges -- dict: id --> Edge

    Vertices and edges are unique and may represent
    several original vertices or edges

    Keep it as simple as possible: first, compute all
    vertices (segment ends and intersections).

    Then, compute the segment-vertex relationships and
    create edges between vertices where one or more
    segments travel.

    The cells are: inside is on the *right* side of
    the half-edges in the cell.

    A cell fully contained in another one has no
    "official" relationship in this class.
    """

    def _add_to_edge(self, v0, v1, sid):
        if (v0, v1) not in self.half_edges:
            self.half_edges[v0, v1] = _HalfEdge((v0, v1))
            self.half_edges[v1, v0] = _HalfEdge((v1, v0))
            self.vertices[v0].add_half_edge((v0, v1))
            self.vertices[v1].add_half_edge((v1, v0))

        self.half_edges[v0, v1].sids.add(sid)

    def __init__(self, sid_segments):  # noqa  XXX needs chopping
        self.half_edges: Dict[Segment, _HalfEdge] = {}
        self.nedges = 0
        self.vertices = {}
        self.cells = {}

        # Map segments out of numpy arrays
        sid_segments = {
            k: (tuple(v[0]), tuple(v[1])) for k, v in sid_segments.items()
        }

        # Remove dups here...
        segments = set()
        segment_sids = collections.defaultdict(list)
        for sid, seg in sid_segments.items():
            simple = tuple(sorted(seg))
            assert simple[0][0].__class__ == FractionClass, (
                simple[0],
                simple[0].__class__,
            )
            assert simple[0][1].__class__ == FractionClass, (
                simple[0],
                simple[0].__class__,
            )
            assert simple[1][0].__class__ == FractionClass, (
                simple[1],
                simple[1].__class__,
            )
            assert simple[1][1].__class__ == FractionClass, (
                simple[1],
                simple[1].__class__,
            )
            segments.add(simple)
            segment_sids[simple].append((sid, seg))

        segments = list(segments)

        global gsegments
        gsegments = segments

        # Find all vertices
        # print('ss', sid_segments)

        vertices = set()
        segment_vertices = collections.defaultdict(set)

        # print('VS')

        def add_vertex(v, segments):
            # if abs(v[0] - -135) < 0.1 and v[0] < -135:
            #     print('AV', fstr(v, segments, approx=True))
            assert len(v) == 2
            # Canonicalize
            v = (Fraction(v[0]), Fraction(v[1]))
            vertices.add(v)
            for segment in segments:
                lsf = line_segment_point_fraction(segment, v)
                assert lsf is not None and 0 <= lsf <= 1

                segment_vertices[segment].add(v)
            return v

        for i, segment in enumerate(segments):
            # print(segment)
            # print(i)
            assert np.any(segment[0] != segment[1])
            add_vertex(segment[0], [segment])
            add_vertex(segment[1], [segment])

        for (inter_v, inter_segments) in scan2d.all_segment_intersections(
            segments
        ):
            add_vertex(inter_v, inter_segments)

        #  inters = list(all_segment_intersections_fast2(segments))

        #  for i, j, inter in inters:
        #      # print(i, j)
        #      if isinstance(inter, SegmentIntersectionPoint):
        #          # print('PT')
        #          add_vertex(inter.v, [segments[i], segments[j]])
        #      else:
        #          # print('LIN')
        #          assert isinstance(inter, SegmentIntersectionLinear), inter
        #          # print(inter)
        #          add_vertex(inter.vs[0], [segments[i], segments[j]])
        #          add_vertex(inter.vs[1], [segments[i], segments[j]])

        self.aabb = AABB()

        for v in vertices:
            self.vertices[v] = _Vertex(v)
            self.aabb.include_point(v)

        # Find all arrangement edges and incidence with orig edges

        for sid, segment in sid_segments.items():
            k = tuple(sorted(segment))
            seg_vertices = segment_vertices[k]

            if False:
                # Slow
                seg_vertices_slow = list()
                for v in vertices:
                    f = line_segment_point_fraction(segment, v)
                    if f is not None and 0 <= f <= 1:
                        seg_vertices_slow.append(v)

                for v in seg_vertices:
                    lsf = line_segment_point_fraction(segment, v)
                    assert lsf is not None and 0 <= lsf <= 1

                assert list(sorted(seg_vertices)) == list(
                    sorted(seg_vertices_slow)
                ), fstr(
                    segment,
                    "K",
                    k,
                    "VERTS",
                    seg_vertices,
                    seg_vertices_slow,
                    approx=True,
                )

                assert line_segment_point_fraction(segment, v) is not None

            v0 = segment[0]

            # Sort going away from v0
            seg_vertices = sorted(
                seg_vertices,
                key=lambda v: (v[0] - v0[0]) ** 2 + (v[1] - v0[1]) ** 2,
            )

            n = len(seg_vertices)
            assert n >= 2
            for i in range(n - 1):
                v0 = seg_vertices[i]
                v1 = seg_vertices[i + 1]

                self._add_to_edge(v0, v1, sid)

        # Check invariants
        for vertex in self.vertices.values():
            he_0 = vertex.half_edges_by_order[0]
            he_cw = he_0
            he_ccw = he_0
            for i in range(len(vertex.half_edges_by_order)):
                # print(i)
                he_cw = vertex.get_next_cw(he_cw)
                he_ccw = vertex.get_next_ccw(he_ccw)
                if i < len(vertex.half_edges_by_order) - 1:
                    assert he_cw != he_0, (
                        he_cw,
                        he_0,
                        vertex.half_edges_by_order,
                    )
                    assert he_ccw != he_0
            assert he_cw == he_0
            assert he_ccw == he_0

        self.general_vector = self._find_general_vector()
        self.general_vector_normal = NP(self.general_vector)[::-1] * (-1, 1)

        self.qtree = None
        self.__init__gstree()

        self.__init__find_cells()

    def __init__find_cells(self):

        # Find all cells
        unprocessed_half_edges = set(self.half_edges.keys())

        # print(fstr(unprocessed_half_edges))

        while len(unprocessed_half_edges) > 0:
            # hyp.note(fstr('UP', unprocessed_half_edges))
            v0, v1 = unprocessed_half_edges.pop()
            # Find the next clockwise edge on v1
            vertices = [v0, v1]

            prev = v0
            cur = v1
            while True:
                prev, cur = vertices[-2:]
                # hyp.note(fstr('V', self.vertices[cur].half_edges_by_order))
                nxt = self.vertices[cur].get_next_ccw((cur, prev))
                # hyp.note(fstr('PCN', prev, cur, nxt))
                if nxt[0] == v0 and nxt[1] == v1:
                    # Complete
                    assert vertices.pop() == v0
                    break
                assert nxt in unprocessed_half_edges
                unprocessed_half_edges.remove(nxt)
                assert nxt[0] == cur
                assert nxt[1] != cur
                assert nxt[1] != prev, fstr(
                    vertices, nxt, self.vertices[cur].half_edges_by_order
                )
                # E.g. for a bowtie, the outer cell
                # will contain the middle vertex twice
                # so cannot assert uniqueness here
                vertices.append(nxt[1])

            cell = tuple(vertices)
            # print(fstr('CELL', cell))
            cell_obj = _Cell()
            area = simple_polygon_area(cell)
            assert area != 0
            cell_obj.negative = area < 0
            cell_obj.count_from_outside = None
            cell_obj.container = None
            cell_obj.containeds = []
            cell_obj.vertices = cell
            cell_obj.contains = []
            self.cells[cell] = cell_obj

            n = len(cell)
            for i in range(n):
                he = self.half_edges[cell[i], cell[(i + 1) % n]]
                assert he.cell is None
                he.cell = cell

        self.__init__find_cell_points()
        # self.__init__find_cell_hierarchy() # XXX NEW
        # self.__init__find_cell_points_new() # XXX # XXX NEW

        # Add cell points to AABB to avoid errors
        # when assuming things with it
        for cell in self.cells.values():
            self.aabb.include_point(cell.point_in_cell)

        self.__init__find_negative_cells_outside()

    def __init__find_cell_hierarchy(self):  # noqa XXX
        """Set cell.contains and cell.container.

        For positive cells, just eat neighbours until
        we find a negative cell with the same half-edge.

        For negative cells, cast a ray.
        """
        done = set()

        print("FCH")
        for cell in self.cells.values():
            print(fstr(cell.vertices, approx=True))
            if cell in done:
                continue
            # Cast ray
            if cell.negative:
                pt = self._find_point_in_cell(cell.vertices)
                loc = locate_point_polygon_winding_number(
                    cell.vertices, cell.point_in_cell
                )
                assert loc < 0
                inters = self.ray_halfedge_query_the_general_ray(
                    pt, self.general_vector
                )
                cells = {}

                for (t, half_edge) in inters:
                    if t < 0:
                        continue
                    assert t != 0
                    # print('HECELL', self.half_edges[half_edge].cell)
                    other_cell = self.half_edges[half_edge].cell
                    if other_cell in cells:
                        del cells[other_cell]
                    else:
                        cells[other_cell] = t

                # Choose the positive ones

                cells = {
                    other_cell: t
                    for other_cell, t in cells.items()
                    if not self.cells[other_cell].negative
                }
                if len(cells) == 0:
                    cell.container = None
                else:
                    cont = list(
                        sorted(cells.items(), key=lambda item: item[1])
                    )[0]
                    cell.container = self.cells[cont[0]]
            else:
                # Positive
                negatives = set()
                positives = {cell}
                doing = list(edges_iter(cell.vertices))
                while doing:
                    v0, v1 = doing.pop()
                    other_vertices = self.half_edges[v1, v0].cell
                    other = self.cells[other_vertices]

                    if other.negative:
                        negatives.add(other)
                    else:
                        if other in done:
                            continue
                        if other not in positives:
                            done.add(other)
                            for half_edge in edges_iter(other.vertices):
                                doing.append(half_edge)
                assert len(negatives) == 1, fstr(
                    (negatives, positives), approx=True
                )
                negative = list(negatives)[0]
                for pos in positives:
                    pos.container = negative

            done.add(cell)

        # Figure out the nesting level of each cell.
        def recurse(cell):
            if cell.count_from_outside is not None:
                return
            if cell.container is None:
                assert cell.negative
                cell.count_from_outside = 0
            else:
                recurse(cell.container)
                cell.count_from_outside = cell.container.count_from_outside + 1

        for cell in self.cells.values():
            recurse(cell)
            if cell.container is not None:
                cell.container.containeds.append(cell)

    def __init__find_cell_points_new(self):
        for cell in self.cells.values():
            if not cell.negative:
                all_points = list(cell.vertices)
                for cont in cell.containeds:
                    all_points.extend(cont.vertices)
                all_points = list(sorted(all_points))
                all_points[0].x

        for cell in self.cells.values():
            if cell.negative:
                if cell.container is None:
                    mp = min(cell.vertices)
                    cell.point_in_cell = (mp[0] - 1, mp[1])
                else:
                    cell.point_in_cell = cell.container.point_in_cell

    def __init__find_cell_points(self):
        """Return a dict from cell to a point in it.

        Useful as a separate function because might
        have optimizations
        """
        for cell in self.cells.values():
            cell.point_in_cell = self._find_point_in_cell(cell.vertices)
            loc = locate_point_polygon_winding_number(
                cell.vertices, cell.point_in_cell
            )
            if cell.negative:
                assert loc < 0
            else:
                assert loc > 0

    def __init__find_negative_cells_outside(self):
        containers = {}
        degrees = {}
        for cell in self.cells.values():
            pt = cell.point_in_cell

            cells = set()
            inters = self.ray_halfedge_query_the_general_ray(
                pt, self.general_vector
            )
            # print('INTERS', inters)
            for (t, half_edge) in inters:
                if t < 0:
                    continue
                assert t != 0
                # print('HECELL', self.half_edges[half_edge].cell)
                cells ^= {self.half_edges[half_edge].cell}

            if not cell.negative:
                # Positive should always contain itself
                assert cell.vertices in cells, (
                    fastr(cell.vertices, cells),
                    cell.vertices,
                    inters,
                    pt,
                )
                cells.remove(cell.vertices)
            else:
                assert cell.vertices not in cells, fstr(
                    cell.vertices, cells, approx=True
                )
            # But it doesn't count
            containers[cell] = list(cells)

            pos_cells = [c for c in cells if not self.cells[c].negative]

            degree = 2 * len(pos_cells) + (0 if cell.negative else 1)
            degrees[cell.vertices] = degree
            cell.count_from_outside = degree

        for cell in self.cells.values():
            cont = containers[cell]
            if len(cont) == 0:
                cell.container = None
                continue
            deg = [degrees[c] for c in cont]
            closest = cont[np.argmax(deg)]
            assert degrees[closest] == degrees[cell.vertices] - 1, fastr(
                degrees[closest],
                degrees[cell.vertices],
                closest,
                cell.vertices,
            )
            assert self.cells[closest].negative != cell.negative
            cell.container = closest

    def __init__qtree(self):
        self.qtree = pmr_quadtree.SegmentPMRQuadTree(
            items=self.half_edges.keys(), threshold=12
        )

    def __init__gstree(self):
        gsegs = []
        for half_edge in self.half_edges.keys():
            # Only do each one once
            if half_edge[0] < half_edge[1]:
                p0 = np.dot(self.general_vector_normal, half_edge[0])
                p1 = np.dot(self.general_vector_normal, half_edge[1])
                seg1d = tuple(sorted((p0, p1)))
                gsegs.append((seg1d[0], seg1d[1], half_edge))

            else:
                assert half_edge[0] > half_edge[1]

        # Stores each half-edge only once
        self.gstree = segment_tree.SegmentTree1D(gsegs)

    def _find_general_vector(self):
        """Generate a (non-unit) vector that is general relative to edges"""

        def key(v):
            return vector_frangle(v)

        frangles = SortedList()

        for half_edge in self.half_edges.values():
            frangle = half_edge.frangle
            frangles.add((frangle + 0 * FRANGLE_90) % MAX_FRANGLE)
            frangles.add((frangle + 1 * FRANGLE_90) % MAX_FRANGLE)
            frangles.add((frangle + 2 * FRANGLE_90) % MAX_FRANGLE)
            frangles.add((frangle + 3 * FRANGLE_90) % MAX_FRANGLE)

        if len(frangles) == 0:
            return (Fraction(0), Fraction(1))

        first = list(islice(frangles.irange(None, None), 1))
        assert len(first) == 1
        first = first[0]

        second = list(
            islice(frangles.irange(first, None, inclusive=(False, True)), 1)
        )
        assert len(second) == 1
        second = second[0]

        assert first != second

        return frangle_unit_square_vector((first + second) / 2)

    def ray_halfedge_query_the_general_ray(
        self, pt, direction, keep_all=False
    ):
        frangle = vector_frangle(direction)
        inters = []

        # Take it outside
        const = 10 * (
            sum(abs(self.aabb.upper - self.aabb.lower))
            / max(abs(NP(direction)))
        )

        rayseg = (pt, NP(pt) + const * NP(direction))
        assert not self.aabb.contains(rayseg[1]), fstr(
            (rayseg, self.aabb), approx=True
        )

        # Special data structure

        offset = np.dot(self.general_vector_normal, pt)

        for low, high, he_key_sorted in self.gstree.stab(offset):
            for he_key in (
                he_key_sorted,
                (he_key_sorted[1], he_key_sorted[0]),
            ):
                half_edge = self.half_edges[he_key]

                self._ray_halfedge_test(
                    half_edge, rayseg, frangle, inters, keep_all, const
                )

        return list(sorted(inters))

    def ray_halfedge_query_general_ray(self, pt, direction, keep_all=False):
        if np.all(direction == self.general_vector) or np.all(
            -direction == self.general_vector
        ):
            return self.ray_halfedge_query_the_general_ray(
                pt, direction, keep_all
            )

        if self.qtree is None:
            raise Exception("Do not want qtree")
            self.__init__qtree()

        frangle = vector_frangle(direction)
        inters = []

        # Take it outside
        const = sum(abs(self.aabb.upper - self.aabb.lower)) / max(
            abs(NP(direction))
        )

        rayseg = (pt, NP(pt) + const * NP(direction))
        assert not self.aabb.contains(rayseg[1]), fstr(
            (rayseg, self.aabb), approx=True
        )

        for he_key in self.qtree.find(rayseg):
            half_edge = self.half_edges[he_key]

            self._ray_halfedge_test(
                half_edge, rayseg, frangle, inters, keep_all, const
            )

        return list(sorted(inters))

    def ray_halfedge_query_general_ray_slow(
        self, pt, direction, keep_all=False
    ):
        """Return all half-edges the given ray intersects.

        The ray must not be parallel to any edge.

        Inclusion of the endpoints:

        If half-edge is cw from direction, startpoint counts,
        else endpoint counts.

        Returns: list of tuples (t, halfedge)
        where t is the solution of pt + t * direction is on
        halfedge.

        Includes negative t.
        """
        frangle = vector_frangle(direction)

        rayseg = (pt, np.array(pt, object) + direction)

        inters = []

        for half_edge in self.half_edges.values():
            self._ray_halfedge_test(
                half_edge, rayseg, frangle, inters, keep_all
            )

        return list(sorted(inters))

    def _ray_halfedge_test(
        self, half_edge, rayseg, frangle, inters, keep_all, const=1
    ):

        (pt, the, tray) = segment_segment_general_intersection(
            half_edge.vs, rayseg
        )

        keep = False

        if keep_all:
            keep = 0 <= the <= 1
        else:
            if 0 < the < 1:
                keep = True
            else:
                # Check direction
                dfrangle = (half_edge.frangle - frangle) % MAX_FRANGLE
                assert dfrangle != 0 and dfrangle != MAX_FRANGLE / 2, (
                    frangle,
                    half_edge.frangle,
                    dfrangle,
                )
                if dfrangle < MAX_FRANGLE / 2:
                    include_start = True
                else:
                    include_start = False

                if the == 0:
                    keep = include_start
                elif the == 1:
                    keep = not include_start

        if keep:
            inters.append((tray * const, half_edge.vs))

    def cell_topology(self):
        """Return a tree of which cell is inside which.

        [cell, [inside_cell, [...]], [inside_cell, ..]]
        """

    def _find_point_in_cell(self, cell):
        """
        ...complex because a cell could be totally contained
        in another.
        """
        half_edge = np.array((cell[0], cell[1]), object)

        pt = (half_edge[0] + half_edge[1]) / 2

        normal = (half_edge[1] - half_edge[0])[::-1] * (1, -1)

        v = self.general_vector
        d = np.dot(v, normal)
        assert d != 0
        if d < 0:
            v = (-v[0], -v[1])
        v = np.array(v, object)

        inters = (
            self.ray_halfedge_query_general_ray
            # self.ray_halfedge_query_general_ray_slow
            (pt, v, keep_all=True)
        )

        # print(fstr(inters, approx=True))

        pos_inters = list(islice(((t, he) for t, he in inters if t > 0), 1))

        if len(pos_inters) == 0:
            return pt + v

        t, he = pos_inters[0]

        global fpc
        fpc = locals()

        return pt + t / 2 * v

    def cell_points(self):
        """Return a dict from cell to a point in it.

        Useful as a separate function because might
        have optimizations
        """
        return {
            cell.vertices: cell.point_in_cell for cell in self.cells.values()
        }

    def poly_cells(self, cells):  # noqa XXX
        """Get the polygons corresponding to the union of the given cells.

        The cells need to be chosen a rule:
        - negative cells can only be included if they
          are canceled by positive cells around them.
        This is how holes are handled.
        """
        global plocs
        plocs = locals()
        # print('POLYCELLS', len(cells), fstr(cells, approx=True))
        m = cell_merger.CellMergerSlow()
        for cell in cells:
            m.add_cell(cell)

        # m.check()

        cells_originals = m.get_cells_originals()
        merged_cells = list(cells_originals.keys())

        positives = []
        negatives = []
        for merged_cell in merged_cells:
            a = simple_polygon_area(merged_cell)
            assert a != 0
            if a > 0:
                positives.append(merged_cell)
            else:
                negatives.append(merged_cell)

        # print('PNCO', fstr(positives, negatives,
        #   cells_originals, approx=True))

        pos_original_cell = {}
        for pos_cell in positives:
            for original in cells_originals[pos_cell]:
                assert original not in pos_original_cell
                pos_original_cell[original] = pos_cell
                # print('POSOR', fastr(original, pos_cell))

        # XXX
        # neg_original_cell = {}
        # for neg_cell in negatives:
        #     for original in cells_originals[neg_cell]:
        #         assert original not in neg_original_cell
        #         neg_original_cell[original] = neg_cell
        #         # print('NEGOR', fastr(original, neg_cell))
        #
        # for cell in cells:
        #     assert (cell in pos_original_cell
        #               or cell in neg_original_cell), fastr(cell)
        #     # XXX Annihilated

        pos_negatives = collections.defaultdict(set)
        if False:
            # This is not correct.
            for negative in negatives:
                originals = list(cells_originals[negative])
                counts = [
                    self.cells[original].count_from_outside
                    for original in originals
                ]
                outermost = originals[np.argmin(counts)]
                cell = self.cells[outermost]
                prevs = set()
                container = cell.container
                assert container is not None
                while (
                    container is not None
                    and container not in pos_original_cell
                ):
                    container = self.cells[container].container
                    assert container not in prevs, prevs
                    prevs.add(container)
                if container is None:
                    # Negative cells can also be not contained within a cell
                    # Find via brute force
                    continue
                outer = pos_original_cell[container]
                pos_negatives[outer].add(negative)
        else:
            for negative in negatives:
                # print('NEG', fstr(negative,
                #    cells_originals[negative], approx=True))
                # All the originals should be in the same outer cell
                outer = None
                for orig in cells_originals[negative]:
                    opos = pos_original_cell.get(orig, None)
                    # print('OPOS', fstr(opos, approx=True))
                    # For negative cells, get the container
                    if opos is None:
                        # Traverse outwards
                        cur = orig
                        while True:
                            cur = self.cells[cur].container
                            # print('Trav', fstr(cur, approx=True))
                            if cur is None:
                                # This happens if there is an outer
                                # ring in two parts
                                break
                            opos = pos_original_cell.get(cur, None)
                            # print('TravO', fastr(opos))
                            if opos is not None:
                                break

                        # ocell = self.cells[cur]
                        # if ocell.negative:
                        #    ocell = ocell.container
                        #    opos = pos_original_cell.get(ocell, None)
                    if opos is not None:
                        if outer is not None:
                            assert outer == opos
                        else:
                            outer = opos

                global locs
                locs = locals()
                assert outer is not None
                pos_negatives[outer].add(negative)

                # pt = (NP(negative[0]) + negative[1]) / 2
                # outer = None
                # for positive in positives:
                #    if locate_point_convex_polygon(00jjkk

        res = []
        for positive in positives:
            # Reverse to be cw
            res.append(
                geom2d_simple_holes.Simple_Polygon_With_Holes(
                    positive,
                    [
                        list(reversed(negverts))
                        for negverts in pos_negatives[positive]
                    ],
                )
            )
        return res
        # XXX

        # Determine containment of cells
        v = self.general_vector

        # hyp.note(fstr('GENV', v))

        cells_for_he = collections.defaultdict(lambda: set())

        for cell, originals in cells_originals.items():
            assert len(cell) >= 3
            m = self.cells[originals[0]].point_in_cell
            #  self.find_point_in_cell(originals[0])

            hes = self.ray_halfedge_query_general_ray(m, v)
            # hyp.note(fstr('QUE', m, cell, 'ORIG', originals))
            for t, he in hes:
                # hyp.note(fstr('HE', t, he))
                assert t != 0
                if t >= 0:
                    assert cell not in cells_for_he[he]
                    cells_for_he[he].add(cell)

        # print(cells_for_he)

        counts = collections.defaultdict(
            lambda: collections.defaultdict(lambda: 0)
        )

        for cell in merged_cells:
            for he in edges_iter(cell):
                for other_cell in cells_for_he[he]:
                    # hyp.note(fstr('Contrib', he, cell, other_cell))
                    if other_cell != cell:
                        counts[other_cell][cell] += 1

        # hyp.note(fstr('COUNTS', counts))

        inside = {}

        for cell in merged_cells:
            inside[cell] = set()
            for possibly_containing_cell, count in counts[cell].items():
                if count % 2 == 1:
                    assert cell != possibly_containing_cell
                    inside[cell].add(possibly_containing_cell)

        # print('inside', inside)

        # Return simple polygons with holes.

        done = set()

        polys = {}

        for k, v in sorted(inside.items(), key=lambda it: len(it[1])):
            # hyp.note(fstr('POLYS', k, v))
            relevant_containers = v - done
            if len(relevant_containers) == 0:
                polys[k] = set()
            else:
                assert len(relevant_containers) == 1
                polys[list(relevant_containers)[0]].add(k)
            done.add(k)

        # print(polys)
        # print('POLYCELLS OUT')

        return [
            geom2d_simple_holes.Simple_Polygon_With_Holes(k, v)
            for k, v in polys.items()
        ]
