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

"""Create a sheetplex from a 3D model
"""

import argparse
import shapely
import shapely.ops
import functools
import operator
import pyrsistent as pyr
import numpy as np
import logging
import warnings

from .sheetplex import SheetPlex, Sheet, Intersection, InterSide
from .plytypes import Fraction
from .geometry.geom2d import Geom2D

logger = logging.getLogger(__name__)


def create_sheetplex(tmesh, params):
    """Create a SheetPlex out of a trimesh mesh.
    """
    # Only stage where we use floats...
    thickness = float(params['thickness'])

    # Turn off Visibledeprecationwarning
    # from trimesh <-> numpy interaction
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',
                                category=np.VisibleDeprecationWarning)

        sheets = create_sheets(tmesh, thickness,
                               offset_tolerance=0.1,
                               normal_tolerance=0.001,
                               only=params['only_sheets'])

    logger.info(f"Found {len(sheets)} sheets")

    sp = SheetPlex(sheets=pyr.m(**sheets), intersections=pyr.m())

    sp = create_intersections(sp)

    return sp


def create_sheets(mesh,  # noqa XXX Too complex function
                  thickness,
                  offset_tolerance,
                  normal_tolerance,
                  only=None):
    """
    Find the sheets that the mesh consists of.

    Heuristically finds the laser-cut sheets needed to construct
    the given mesh.

    Returns
        dict of sheet id -> Sheet object
    """
    if only is not None:
        only = set(only)

    sheets = {}

    idx = 0
    all_facets = find_all_facets(mesh)
    todo_facets = {i for i in np.arange(len(all_facets.facets))}
    facet_order = reversed(np.argsort(all_facets.areas))

    n_slices = 3

    for facet in facet_order:
        if facet in todo_facets:
            todo_facets.remove(facet)

            normal, offset = _facet_plane(all_facets, facet)
            logger.debug(f'Facet {facet}: {normal} {offset}')

            oppnormal, oppoffset = -normal, -(offset - thickness)
            sides = np.array([
                list(normal) + [offset],
                list(oppnormal) + [oppoffset],
            ], np.float64)

            out_tolerance = offset_tolerance
            slices_height = [offset + out_tolerance,
                             offset - thickness - out_tolerance] + list(
                np.linspace(offset - thickness + offset_tolerance,
                            offset - offset_tolerance, n_slices))

            slices = mesh.section_multiplane([0, 0, 0], normal, slices_height)
            outside_slices = [Geom2D.from_shapely(
                shapely.ops.unary_union(sl.polygons_full))
                if sl else Geom2D.empty()
                # sl if sl else Geom2D.empty()
                for sl in slices[0:2]]
            slices = slices[2:]

            for sl in slices:
                assert sl is not None, (
                    "Null slice - wrong thickness setting?",
                    normal, offset)

            slice_rotation = slices[0].metadata['to_3D'].copy()
            # Remove translation
            slice_rotation[0:3, 3] = 0

            slices = [Geom2D.from_shapely(
                shapely.ops.unary_union(sl.polygons_full))
                for sl in slices]

            same = [facet]
            opposed = []
            within = []

            for other_facet in list(todo_facets):
                assert other_facet != facet
                _find_facet_status(other_facet, same, opposed, within,
                                   mesh,
                                   all_facets,
                                   todo_facets,
                                   normal, offset,
                                   thickness,
                                   normal_tolerance, offset_tolerance)
            if len(opposed) == 0:
                logger.debug('Skipping one at %s, no opposed facets', idx)
                continue

            id = 's%d' % idx
            idx += 1

            logger.debug(f'Adding sheet {id}')

            # Find area supported by actual faces

            trb = slice_rotation.T
            faces = []
            for fs in (same, opposed):
                fsupport = Geom2D.empty()
                for f in fs:
                    for tri in all_facets.facets[f]:
                        verts = trb[:3, :3].dot(
                            mesh.triangles[tri].T).T[:, 0:2]
                        verts = vFraction(verts)
                        fsupport = fsupport | Geom2D.polygon(
                            verts, reorient=True)
                faces.append(fsupport)

            slices_max = functools.reduce(operator.or_, slices)
            slices_min = functools.reduce(operator.and_, slices)

            if only is not None and id not in only:
                logger.debug('Skipping due to "only" setting')
                continue

            sheets[id] = Sheet(
                id=id,
                sides=sides,
                slice_rotation=slice_rotation,
                inverse_slice_rotation=slice_rotation.T,
                faces=faces,
                outside_slices=outside_slices,
                both_faces=faces[0] & faces[1],
                slices_max=slices_max,
                slices_min=slices_min
            )

    return sheets


def _find_facet_status(other_facet, same, opposed, within,
                       mesh,
                       all_facets,
                       todo_facets,
                       normal, offset,
                       thickness,
                       normal_tolerance, offset_tolerance):
    other_normal, other_offset = _facet_plane(all_facets, other_facet)

    norm = np.linalg.norm

    # Same plane?
    if norm(other_normal - normal) < normal_tolerance and \
            norm(other_offset - offset) < offset_tolerance:
        same.append(other_facet)
        todo_facets.remove(other_facet)
        return

    # Opposing side of plane
    if norm(-other_normal - normal) < normal_tolerance and \
            norm(-other_offset - offset + thickness) < \
            offset_tolerance:
        opposed.append(other_facet)
        todo_facets.remove(other_facet)
        return

    # Within plane
    if abs(np.dot(other_normal, normal)) < normal_tolerance:
        # Check that all points are within this and other
        boundary_verts = set()
        for v0, v1 in all_facets.boundaries[other_facet]:
            boundary_verts.add(v0)
            boundary_verts.add(v1)
        is_within = True
        for v in boundary_verts:
            if not (offset - thickness - offset_tolerance
                    <= np.dot(mesh.vertices[v], normal)
                    <= offset + offset_tolerance):
                is_within = False
                break
        if is_within:
            within.append(other_facet)
            todo_facets.remove(other_facet)
            return


vFraction = np.vectorize(Fraction, [object])
"""Vectorize a numpy array"""


def _facet_plane(all_facets, facet):
    """Get the plane dot vector and offset for a facet.
    """
    origin = all_facets.origins[facet]
    normal = all_facets.normals[facet]
    return normal, np.dot(origin, normal)
#


def find_all_facets(mesh):
    """Return all facets of a mesh.

    Work around a trimesh API pain point where only
    sides with 2 or more triangles count as facets.

    Returns: a Namespace with the following parallel facet lists
        facets - list of lists of triangle ids
        normals - list of normal vectors
        origins - list of a single point on each facet
        areas - list of total area of the facet
        boundaries - list of facet boundary vertices
    """
    assert len(mesh.faces) == len(mesh.triangles)
    triangles_left = set(range(len(mesh.triangles)))

    facets = list(mesh.facets)
    normals = list(mesh.facets_normal)

    if mesh.facets_origin is not None:
        # For some reason, this gets None and not []
        origins = list(mesh.facets_origin)
    else:
        origins = []

    areas = list(mesh.facets_area)
    boundaries = list(mesh.facets_boundary)
    for facet in facets:
        triangles_left -= set(facet)

    for triangle_id in triangles_left:
        facets.append([triangle_id])
        normals.append(mesh.face_normals[triangle_id])
        areas.append(mesh.area_faces[triangle_id])
        origins.append(mesh.triangles[triangle_id][0])
        b = []
        for i in range(3):
            b.append((mesh.faces[triangle_id][i],
                      mesh.faces[triangle_id][(i + 1) % 3]))
        boundaries.append(b)

    return argparse.Namespace(
        facets=facets,
        normals=normals,
        origins=origins,
        areas=areas,
        boundaries=boundaries,
    )


def create_intersections(sp):
    sheets = list(sp.sheets.values())
    for i, sheet0 in enumerate(sheets):
        for sheet1 in sheets[i+1:]:
            sp = create_intersection(sp, (sheet0, sheet1))
    return sp


def create_intersection(sp, sheets):
    id = intersection_id(sheets)

    sheet1, sheet2 = sheets

    # Intersection line direction in 3-space
    inter_dir = np.cross(sheet1.normal(0), sheet2.normal(0))
    ilen = np.linalg.norm(inter_dir)

    if ilen < 0.0001:  # Parallel -- no intersection
        return sp

    inter_dir /= ilen

    # Intersection origin in 3-space must be on both planes.
    # Lstsq chooses such point closest to 3D origin.
    corner_origs = []
    for t1 in (0, 1):
        for t2 in (0, 1):
            inter_orig, res, rank, s = np.linalg.lstsq(
                [sheet1.normal(t1), sheet2.normal(t2)],
                [sheet1.offset(t1), sheet2.offset(t2)],
                rcond=None
            )
            corner_origs.append(inter_orig)

# XXX Future
#    # Figure out shortening of teeth
#    d0 = np.linalg.norm(corner_origs[0] - corner_origs[3])
#    d1 = np.linalg.norm(corner_origs[1] - corner_origs[2])
#
#    if d0 > params['max_finger_length']:
#    elif d1 > parmas['max_finger_length']:

    inter = Intersection(
        id=id,
        direction=inter_dir,
        origins=corner_origs,
    )

    inter_sides = [create_inter_side(sp, inter, sheet, idx) for
                   (idx, sheet) in enumerate(sheets)]

    inter = inter.set('sides', tuple(inter_sides))

    return sp.transform(('intersections',),
                        lambda inters: inters.set(id, inter))


def intersection_id(sheets):
    return ",".join(list(sorted([sheet.id for sheet in sheets])))


def create_inter_side(sp, inter, sheet, idx):
    id = (inter.id, sheet.id)

    direction = sheet.project_point4(
        list(inter.direction) + [0.]).astype(np.float64)
    direction = vFraction(direction)
    normal = direction[::-1] * [1, -1]

    sheet_inter_origs = [
        sheet.project_point4(list(orig) + [1.])
        for orig in inter.origins
    ]
    sheet_inter_origs = vFraction(sheet_inter_origs)

    normal_offsets = [np.dot(orig, normal) for orig in sheet_inter_origs]

    return InterSide(
        id=id,
        intersection=inter.id,
        sheet=sheet.id,
        idx=idx,
        direction=direction,
        normal=normal,
        origin_offset=direction.dot(sheet_inter_origs[0]),
        min_normal_offset=np.min(normal_offsets),
        max_normal_offset=np.max(normal_offsets),
    )
