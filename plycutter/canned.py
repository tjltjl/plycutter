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

"""
Canned routines to run a bunch of heuristic steps for generating finger joints.
"""

import logging
import argparse

from . import sheetbuild
from . import heuristics as heu

logger = logging.getLogger(__name__)

sb_log = None
"""A global variable where the log from canned_1 is stored for debugging.

Use this variable in jupyterlab to poke around.
"""


def canned_1(sheetplex, params, show=False):
    """Run a number of heuristic steps to makek a finger jointed 2d pattern"""

    global sb_log
    sb_log = []
    ok = False
    try:
        sp = sheetplex
        logger.info("Create sheetbuild")
        sb = sheetbuild.create_sheetbuild(sp, params)
        sb_log.append(sb.set('dbg', 'create'))

        sb.check()

        if show:
            do_show('Creation', sb)

        logger.info("Process")
        for op in [
            heu.choose_obviously_irrelevant,
            heu.remove_loose,
            heu.update_intersection_ok,
            heu.update_sheet_ok,
            heu.update_intersection_ok,
            heu.remove_loose,
            heu.update_intersection_ok,
            heu.heuristic_multi_inter_single_decisions,
            heu.update_intersection_ok,
            heu.choose_the_ok_side,
            heu.single_fingers_random,
            heu.clean_up_thin_ok_parts,
            # select_all,
        ]:
            if op is None:
                break
            logger.info(op)
            sb_log.append(sb.set('dbg', op))
            sb = op(sb, params)
            if show:
                do_show(op, sb)

        sb_log.append(sb)
        sb.check(force=True)

        ok = True

    except Exception as e:
        logger.error('Exception in canned_1: %s', e, exc_info=True)
        if not params['return_on_failure']:
            raise

    return argparse.Namespace(
        ok=ok,
        sheetbuild=sb,
        history=sb_log,
        sheet_cuts=sb.sheet_chosen,
        params=params,
    )


def do_show(stage, sb):
    from IPython.display import display
    from pylab import subplots
    sp = sb.sheetplex
    s = list(sorted(sp.sheets.keys()))
    N = 4
    print('AFTER STAGE', stage)
    i = 0
    fig = None
    for sheet_id in s:
        if i % N == 0:
            if fig is not None:
                display(fig)
            fig, ax = subplots(1, N, figsize=(15, 4))
        ax[i % N].set_title(sheet_id)
        sheetbuild.show_sheet(ax[i % N], sb.sheetplex, sb, sheet_id)
        i += 1
    if fig is not None:
        display(fig)
