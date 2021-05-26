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

from ..geom2d import Geom2D

def test_buffer():
    square = Geom2D.rectangle(0, 0, 1, 1)
    bsquare = square.buffer(1, resolution=16)
    assert bsquare.locate((0.5, 0.5)) == 1
    assert bsquare.locate((1.5, 0.5)) == 1
    assert bsquare.locate((2.5, 0.5)) == -1
