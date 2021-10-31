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

import argparse
import logging
import pathlib
import sys

import trimesh

from . import create_sheetplex
from . import canned
from . import writer
from .plytypes import F

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

logging.getLogger('trimesh').setLevel(logging.INFO)
logging.getLogger('ezdxf').setLevel(logging.WARN)
logging.getLogger('chardet.charsetprober').setLevel(logging.WARN)
logging.getLogger('matplotlib').setLevel(logging.WARN)


def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--thickness', type=F, default=6,
                        help='Set the thickness of sheets to find.')
    parser.add_argument('--min_finger_width', type=F, default=3,
                        help='Set minimum width for generated fingers.')
    parser.add_argument('--max_finger_width', type=F, default=5,
                        help='Set maximum width for generated fingers.')
    parser.add_argument('--support_radius', type=F, default=12,
                        help='Set maximum range for generating material '
                        'on a sheet where '
                        'neither surface is visible')
    parser.add_argument('--debug', action="store_true",
                        help='Turn on debugging.')

    parser.add_argument('--vdebug', action="store_true",
                        help='Turn on visual debugging (try in jupyterlab).')

    parser.add_argument('--final_dilation', type=F, default=F(5, 100),
                        help='Final dilation (laser cutter kerf compensation)')

    parser.add_argument('--random_seed', type=int, default=42,
                        help='Random seed for pseudo-random heuristics')

    # XXX Doesn't convert properly yet
    parser.add_argument('--only_sheets', type=str, default=None,
                        help='Not implemented yet')

    parser.add_argument('-o', '--output_file', type=str, default=None,
                        help='File to write the output in')
    
    parser.add_argument('-f', '--format', type=str, default="dxf",
                        help='dxf or svg output')
    
    parser.add_argument('infile', type=str,
                        help='STL file to process')

    args = parser.parse_args(arguments)

    infile = pathlib.Path(args.infile)

    args.format = args.format.lower()
    if args.format not in ('svg', 'dxf'):
        raise Exception("Improper file format")

    if args.output_file is None:
        outfile = infile.with_suffix("." + args.format)
    else:
        outfile = args.output_file

    logger.info(f'Loading "{infile}. Will write "{outfile}""')
    mesh = trimesh.load(str(infile))
    logger.info(f'Parameters: {args}')

    logger.info('Creating sheetplex')
    sp = create_sheetplex.create_sheetplex(mesh, vars(args))
    if len(sp.sheets) == 0:
        raise Exception("No sheets found")

    logger.info('Running canned heuristic steps')

    result = canned.canned_1(sp, vars(args), args.vdebug)

    logger.info(f'Writing "{outfile}"')
    to_write = result.sheet_cuts

    to_write_dilated = {
        k: v.buffer(args.final_dilation)
        for k, v in to_write.items()
    }
    print(to_write_dilated)

    if args.format == 'dxf':
        writer.write_dxf(str(outfile), to_write_dilated)
    elif args.format == 'svg':
        writer.write_svg(str(outfile), to_write_dilated)

    logger.info('Done!')
