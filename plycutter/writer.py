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

import ezdxf
import numpy as np
import logging

# XXX
from .geometry.aabb import AABB

logger = logging.getLogger(__name__)


def write_svg(filename, geom2ds):
    file = open(filename, "w")
    file.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    file.write('<svg version = "1.1" xmlns="http://www.w3.org/2000/svg">\n')

    # Trivial nesting horizontally
    x = 0
    for name, geom in geom2ds.items():
        print(name)
        if geom.is_empty():
            logger.warn(f"Empty sheet {name}")
            continue

        aabb = AABB()

        for polygon in geom.polygons():
            coords = np.array(polygon.spwhs[0].outer)
            for pt in coords:
                aabb.include_point(pt)

        margin = 3
        x_offset = x - aabb.lower[0] + margin
        y_offset = 0 - aabb.lower[1]
        x += aabb.upper[0] - aabb.lower[0] + margin

        def draw_coords(coords):
            coords = list(coords)
            coords = np.array(coords + coords[0:1])  # Loop
            coords = coords + [x_offset, y_offset]
            coords = coords.astype(np.float64)
            assert np.all(np.isfinite(coords))
            path = "M "
            for coord in coords:
                path += "{:0.6f},{:0.6f} ".format(coord[0], coord[1])
            style = (
                "fill:none;stroke:#000000;"
                "stroke-width:1px;stroke-opacity:1.0"
            )
            file.write(f"""<path d ="{path}" style="{style}" />\n""")

        for polygon in geom.polygons():
            draw_coords(polygon.spwhs[0].outer)
            for hole in polygon.spwhs[0].holes:
                draw_coords(hole)
    file.write("</svg>")
    file.close()


def write_dxf(filename, geom2ds):
    dwg = ezdxf.new("AC1015")
    modelspace = dwg.modelspace()

    # Trivial nesting horizontally
    x = 0
    for name, geom in geom2ds.items():
        print(name)
        if geom.is_empty():
            logger.warn(f"Empty sheet {name}")
            continue

        aabb = AABB()

        for polygon in geom.polygons():
            coords = np.array(polygon.spwhs[0].outer)
            for pt in coords:
                aabb.include_point(pt)

        margin = 3
        x_offset = x - aabb.lower[0] + margin
        y_offset = 0 - aabb.lower[1]
        x += aabb.upper[0] - aabb.lower[0] + margin

        def draw_coords(coords):
            coords = list(coords)
            coords = np.array(coords + coords[0:1])  # Loop
            coords = coords + [x_offset, y_offset]
            coords = coords.astype(np.float64)
            assert np.all(np.isfinite(coords))
            modelspace.add_lwpolyline(coords.tolist())

        for polygon in geom.polygons():
            draw_coords(polygon.spwhs[0].outer)
            for hole in polygon.spwhs[0].holes:
                draw_coords(hole)

    dwg.saveas(filename)
