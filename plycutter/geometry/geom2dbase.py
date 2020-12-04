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

from abc import abstractmethod
import matplotlib
import matplotlib.path
import matplotlib.patches


class Geom2DBase:
    """A set of semi-open convex polygons.

    (called regularized polygons in CGAL)

    Supports the set operations ``&``, ``|``, ``^``, ``-``
    as well as ``~``..

    Like Geom1D, we get a convenient set of properties
    by asserting that the enclosing lines and
    the contents are in different spaces (
    e.g. rationals and rationals + sqrt(2))
    so that unions and differences work in
    the most natural way possible.
    """
    @classmethod
    def polygon(cls, points, reorient=False):
        pass

    @classmethod
    def union(cls, geoms):
        raise Exception()

    @classmethod
    def rectangle(cls, x0, y0, x1, y1):
        return cls.polygon([
            (x0, y0),
            (x0, y1),
            (x1, y1),
            (x1, y0),
        ], reorient=True)

    @abstractmethod
    def __or__(self, other): pass
    @abstractmethod
    def __and__(self, other): pass
    @abstractmethod
    def __sub__(self, other): pass

    @abstractmethod
    def locate(self, pt):
        """Locate point wrt polyset

        1 = in
        -1 = out
        0 = indeterminate, possibly on edge,
            possibly on internal, virtual edge

        """

    @abstractmethod
    def all_segments(self):
        """Iterate over all line segments defining this geom.
        Depending on representation,
        may include internal segments that are not
        actual external boundary.
        """

    @abstractmethod
    def polygons(self):
        """Iterate over all distinct polygons (with holes)
        in this geom.
        """

    @abstractmethod
    def exterior(self):
        """Get the outer boundary as a list of vertices,
        provided this is a single polygon.
        """

    @abstractmethod
    def holes(self):
        """Iterator over the vertex lists of holes,
        provided this is a single polygon.
        """

    def to_mpatch(self, color='red', alpha=1.0, label=None,
                  linewidth=0, ax=None):
        """Convert to matplotlib.patches.PathPatch"""
        codes = []
        verts = []
        # print('MPATCH', self)

        def add_path(new_verts):
            new_verts = list(new_verts)
            assert len(new_verts) > 1
            new_verts = new_verts + new_verts[0:1]
            codes.extend([matplotlib.path.Path.MOVETO] +
                         (len(new_verts) - 1) * [matplotlib.path.Path.LINETO])
            verts.extend(new_verts)

        for poly in self.polygons():
            ext = poly.exterior()
            add_path(ext)

            for hole in poly.holes():
                # print('HOLE', hole)
                add_path(reversed(hole))

        if len(verts):
            path = matplotlib.path.Path(verts, codes)
            patch = matplotlib.patches.PathPatch(path,
                                                 fill=True,
                                                 facecolor=color,
                                                 edgecolor='black',
                                                 label=label,
                                                 alpha=alpha,
                                                 linewidth=linewidth)
        else:
            patch = None

        if ax is not None:
            ax.add_patch(patch)
            ax.autoscale_view()
        return patch

    def show2d(self, ax, color, alpha=1.0, label=None, linewidth=0):
        """Plot this Geom2D onto the given matplotlib Axes object.
        """
        patch = self.to_mpatch(color=color, alpha=alpha, label=label,
                               linewidth=linewidth)
        if patch is not None:
            ax.add_patch(patch)
        ax.set_aspect(1)
        ax.autoscale_view()
