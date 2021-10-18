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

"""Heuristic creation of finger joints.

This module contains heuristic functions that take
as their parameters a ``SheetBuild`` and a params dict
and return a new version of the ``SheetBuild`` with some heuristic
applied..
"""

import pyrsistent as pyr
import logging
import collections
import numpy as np

from .plytypes import F, fastr
from .misc import Autocompletable
from .sheetbuild import map_and, map_or
from .geometry.geom1d import Geom1D
from .geometry.geom2d import Geom2D

logger = logging.getLogger(__name__)


def choose_obviously_irrelevant(sb, params=None):
    """For joints where there is nothing next to the joint
    on a sheet, unchoose the sheet.
    """
    sb.check()
    sp = sb.sheetplex
    unchoices = {}
    for interside in sp.intersides():
        ok = sb.interside_ok[interside.id]
        fing = interside.make_fingers(ok.buffer(F(1, 100)))
        outside = sp.sheets[interside.sheet].slices_max & fing.buffer(F(1, 100))
        outside -= fing.buffer(F(1, 1000))
        outsideless = ~interside.project_to_1d(outside)
        opposite = interside.opposite(sp)
        unchoices[interside.id] = outsideless & ok & sb.interside_ok[opposite.id]
        # logger.debug(interside.id, outsideless, ok, outsideless & ok)

    # logger.debug(choices)
    return sb.unchoose_intersides(unchoices)


def remove_loose(sb, params=None):
    """Remove loose pieces that do not connect with the
    both-face support"""
    sb.check()
    sp = sb.sheetplex
    for sheet in sp.sheets.values():
        full_support = sheet.faces[0] & sheet.faces[1]
        possible = sb.sheet_chosen[sheet.id] | sb.sheet_ok[sheet.id]
        impossible = Geom2D.empty()
        for polygon in possible.polygons():
            if (polygon & full_support).is_empty():
                # assert (sb.sheet_chosen[sheet.id] & polygon).is_empty()
                assert (sb.sheet_chosen[sheet.id] & polygon).area() < 0.01, (
                    sheet.id,
                    (sb.sheet_chosen[sheet.id] & polygon).area(),
                )
                impossible |= polygon
        sb = sb.transform(["sheet_ok", sheet.id], lambda old: old - impossible)
    return sb


def update_intersection_ok(sb, params=None):
    """A position cannot be ok for an intersection is it's not
    possible to put something there"""
    sb.check()
    for interside in sb.sheetplex.intersides():
        ok = sb.interside_ok[interside.id]
        fing = interside.make_fingers(ok.buffer(F(1, 10)))
        max_fing = fing & sb.sheet_ok[interside.sheet]  # .buffer(-.001)
        unsupported = interside.project_to_1d(fing) - interside.project_to_1d(max_fing)
        sb = sb.transform(
            ["interside_ok", interside.id],
            lambda old: old - unsupported.buffer(F(1, 1000)),
        )
    return sb


def update_sheet_ok(sb, params=None):
    """A point on the sheet cannot be ok if it is on intersections
    but none of the intersides are ok for it.
    """
    sb.check()
    sp = sb.sheetplex
    for sheet in sp.sheets.values():
        on_interside = Geom2D.empty()
        ok_on_interside = Geom2D.empty()
        for interside in sp.intersides(sheet.id):
            opposite = interside.opposite(sp)
            ok_on_interside |= interside.make_fingers(sb.interside_ok[interside.id])
            on_interside |= interside.make_fingers(
                sb.interside_ok[interside.id] | sb.interside_ok[opposite.id]
            )

        sb = sb.unchoose_sheet_area(
            sheet.id,
            ((on_interside - ok_on_interside).buffer(F(2, 1000)))
            - sb.sheet_chosen[sheet.id],
        )
    return sb


class _Interaction(pyr.PRecord, Autocompletable):
    interside = pyr.field()
    range = pyr.field()
    other_interside = pyr.field()
    other_range = pyr.field()


class _Interactions(pyr.PRecord, Autocompletable):
    sheetplex = pyr.field()
    by_intersection = pyr.field()


def _find_interactions(sb):
    """Find the intersections between intersections.

    We want the structure of which parts are interdependent.

    The basis is: region X of intersection A interferes with
    region Y of intersection B.
    """
    sp = sb.sheetplex

    # Set of pairs of intersections that intersect
    # Map: intersection -> [range, sheet, other_interside, other_range]

    interactions = collections.defaultdict(lambda: [])

    for sheet in sp.sheets.values():
        intersides = list(sp.intersides(sheet.id))
        for i, i_interside in enumerate(intersides):
            i_fingers = i_interside.make_fingers(sb.interside_ok[i_interside.id])

            for j_interside in intersides[i + 1 :]:
                j_fingers = j_interside.make_fingers(sb.interside_ok[j_interside.id])
                double = i_fingers & j_fingers
                if not double.is_empty():
                    i_double = i_interside.project_to_1d(double)
                    j_double = j_interside.project_to_1d(double)

                    interactions[i_interside.intersection].append(
                        _Interaction(
                            interside=i_interside.id,
                            range=i_double,
                            other_interside=j_interside.id,
                            other_range=j_double,
                        )
                    )

                    interactions[j_interside.intersection].append(
                        _Interaction(
                            interside=j_interside.id,
                            range=j_double,
                            other_interside=i_interside.id,
                            other_range=i_double,
                        )
                    )

    return _Interactions(
        sheetplex=sb.sheetplex,
        by_intersection=interactions,
    )


def _find_interaction_clusters(sb, interactions):
    sp = sb.sheetplex

    cluster_members = {}

    n = 0

    for intersection_id, interaction_list in interactions.by_intersection.items():
        # logger.debug("\n", intersection_id)
        for interaction in interaction_list:
            cmap = pyr.pmap(
                {
                    sp.interside(interaction.interside).intersection: interaction.range,
                    sp.interside(
                        interaction.other_interside
                    ).intersection: interaction.other_range,
                }
            )

            # Check which clusters it belongs to
            belongs_to = []
            for k, v in cluster_members.items():
                if len(map_and(v, cmap)) > 0:
                    belongs_to.append(k)
            # logger.debug('belongs_to', belongs_to)

            if len(belongs_to) > 0:
                # logger.debug('add_to', belongs_to[0])
                cluster_members[belongs_to[0]] = map_or(
                    cluster_members[belongs_to[0]], cmap
                )

                for k in belongs_to[1:]:
                    cluster_members[belongs_to[0]] = map_or(
                        cluster_members[belongs_to[0]], cluster_members[k]
                    )
                    del cluster_members[k]
                # logger.debug('-->', cluster_members[belongs_to[0]])
            else:
                # If none, start a new cluster
                cid = "c%d" % n
                n += 1
                cluster_members[cid] = cmap
                # logger.debug('New', cid, cluster_members[cid])

    return cluster_members


# noqa because this function is too complex; needs to be redone anyway
# because the heuristic of just choosing a single sheet for a cluster
# is not the best way.
def heuristic_multi_inter_single_decisions(sb, params=None):  # noqa
    """Set each multi-intersection to *one* sheet
    for all intersections, plus
    add enough extra to make it fit.
    """
    logger.debug("Find interactions")
    interactions = _find_interactions(sb)
    logger.debug("Find clusters")
    clusters = _find_interaction_clusters(sb, interactions)

    # Reserve all multi-intersections from
    # being used as the extra connector
    intersection_reserved = pyr.pmap({})
    for cluster in clusters.values():
        intersection_reserved = map_or(intersection_reserved, cluster)

    sp = sb.sheetplex
    for cluster in clusters.values():
        logger.debug("CLUSTER", cluster)
        intersections = [
            sp.intersections[intersection_id] for intersection_id in cluster.keys()
        ]

        sheet_ids = set()
        intersides_by_sheet = collections.defaultdict(lambda: [])

        for intersection in intersections:
            for interside in intersection.sides:
                sheet_ids.add(interside.sheet)
                intersides_by_sheet[interside.sheet].append(interside)

        sheets = [sp.sheets[sheet_id] for sheet_id in sheet_ids]

        eps = F(2, 10)
        min_finger = 3

        # Compute total inserted area for each sheet
        areas = []
        needs_extra = []
        for sheet in sheets:
            fingers = Geom2D.empty()
            for interside in intersides_by_sheet[sheet.id]:
                fingers |= interside.make_fingers(
                    cluster[interside.intersection] & sb.interside_ok[interside.id]
                )
            areas.append(fingers.area())

            # Check whether we need extra connection to the
            # main body. E.g. triple 90-degree joints require this
            # because the corner piece is not automatically connected
            # to the joint body, whereas some other joints may not
            area = (fingers.buffer(eps) & sb.sheet_chosen[sheet.id]).area()
            needs_extra.append(area <= min_finger * eps)

        # logger.debug('CHOOSE', np.argmax(areas),
        #       [sheet.id for sheet in sheets],
        #        areas, needs_extra)

        choice = np.argmax(areas)
        sheet = sheets[choice]

        if needs_extra[choice]:
            found = False
            logger.debug(f"NEEDS EXTRA {choice} {sheet.id} {cluster}")
            for interside in intersides_by_sheet[sheet.id]:
                # logger.debug(interside.id, sb.interside_ok[interside.id],
                #    interside)

                # Can we add it here?
                intersection = sp.intersections[interside.intersection]
                default_range = cluster[intersection.id]
                max_range = default_range.buffer(min_finger + eps)

                # Potential trouble is this is in the middle?
                max_range &= sb.interside_ok[interside.id]
                max_range -= default_range

                extra_ints = sorted(
                    max_range.intervals, key=lambda iv: (iv[1] - iv[0], iv[0])
                )
                logger.debug(f"  EXTRA INTS f{interside.id} f{max_range}")
                # logger.debug(extra_ints)
                if len(extra_ints) == 0:
                    continue

                extra_interval = extra_ints[-1]

                if extra_interval[1] - extra_interval[0] < 0.5 * min_finger:
                    logger.debug("  TOO SMALL")
                    continue

                found = True
                # logger.debug(extra_interval)
                cluster = cluster.transform(
                    [interside.intersection],
                    lambda prev: prev
                    | Geom1D(
                        [
                            extra_interval,
                        ]
                    ),
                )
                break

            if not found:
                raise Exception("NOT FOUND %s" % cluster)
        else:
            logger.debug(f"NO EXTRA f{choice} {sheet.id} f{cluster}")

        # Choose the side of the joint
        # on all intersides on the chosen sheet

        sb.check()
        # logger.debug("   OK %s\n     CHOSEN %s",
        #    sb.interside_ok, sb.interside_chosen)
        logger.debug("Make choices")
        for interside in intersides_by_sheet[sheets[choice].id]:
            logger.debug(
                "Choose: %s %s %s",
                interside.id,
                cluster[interside.intersection],
                cluster,
            )
            chosen_range = (
                cluster[interside.intersection] & sb.interside_ok[interside.id]
            )
            sb = sb.choose_interside(interside.id, chosen_range)

        # Unchoose all other intersides

        for other_intersection_id, rnge in cluster.items():
            other_intersection = sp.intersections[other_intersection_id]
            for side in other_intersection.sides:
                if side.sheet != sheets[choice].id:
                    logger.debug("  Unchoose %s %s", side.id, rnge)
                    sb = sb.unchoose_interside(side.id, rnge)

        logger.debug("   OK %s\n     CHOSEN %s", sb.interside_ok, sb.interside_chosen)
        sb.check()

    return sb


# def resolve_multi_intersections_1(sb, multi):


def choose_the_ok_side(sb, params=None):
    """If only one side of an intersection is ok, choose it.

    To be run only after multi-intersections handled.
    """
    sb.check()

    sp = sb.sheetplex
    for intersection in sp.intersections.values():
        for interside in intersection.sides:
            this_only = sb.interside_ok[interside.id]
            this_only -= sb.interside_chosen[interside.id]
            this_only -= sb.interside_ok[interside.opposite(sp).id]
            this_only -= sb.interside_chosen[interside.opposite(sp).id]

            buf = F(1, 100)
            orig_this_only = this_only
            this_only = this_only.buffer(-buf).buffer(buf)
            small_bits = orig_this_only - this_only

            logger.debug("Choose onlyok: %s %s", interside.id, this_only)
            sb = sb.choose_interside(interside.id, this_only)
            sb = sb.unchoose_interside(interside.id, small_bits)

            sb.check()

    return sb


def random_min_max(n, mn, mx, total, rng):
    assert mn * n <= total <= mx * n
    lengths = np.zeros(n, dtype=np.float64)
    order = np.arange(n)
    rng.shuffle(order)
    running_total = 0
    n_handled = 0
    for i in order:
        left = total - running_total
        # If all the rest are at mx, we can do n
        at_most = left - (n - n_handled - 1) * mn
        at_min = left - (n - n_handled - 1) * mx
        # logger.debug(i, left, at_min)
        low = max(at_min, mn)
        hi = min(at_most, mx)
        # logger.debug(low, hi)
        assert low <= hi
        cur = rng.uniform(low, hi)
        # logger.debug(cur)
        lengths[i] = cur
        running_total += cur
        n_handled += 1
    # logger.debug(lengths)
    result = []
    cur_sum = 0
    for i in range(n):
        result.append((cur_sum, cur_sum + lengths[i]))
        cur_sum += lengths[i]
    assert cur_sum - total < 1e-8, (lengths, sum(lengths), result, cur_sum, total)
    return result


def single_fingers_random(sb, params=None):
    min_width = params["min_finger_width"]
    max_width = params["max_finger_width"]
    seed = params["random_seed"]
    sb.check()
    rng = np.random.RandomState(seed)
    # XXX Assert that there are no multiple intersections
    # any more
    sp = sb.sheetplex
    for intersection in sp.intersections.values():
        sides = intersection.sides
        assert len(sides) == 2
        # At this point, everything left should be
        # ok on both sides
        ok = [sb.interside_ok[side.id] for side in sides]
        chosen = [sb.interside_chosen[side.id] for side in sides]
        available = [ok[i] - chosen[i] for i in range(len(ok))]
        assert (available[0] ^ available[1]).measure1d() < 1e-8, (
            intersection,
            ok[0],
            ok[1],
            available[0],
            available[1],
        )
        # Start assigning by intervals
        for (start, end) in available[0].intervals:
            logger.debug(f"{intersection.id} : {fastr(start, end)}")
            geom1d = Geom1D([(start, end)])
            end_dist = 2
            on_side = []
            for end_interval in [(start - end_dist, start), (end, end + end_dist)]:
                ei = Geom1D([end_interval])
                on_sides = [
                    (sb.interside_chosen[side.id] & ei).measure1d() for side in sides
                ]
                on_side.append(np.argmax(on_sides))

            odd = on_side[0] != on_side[1]

            n_min_base = np.ceil((end - start) / max_width)
            n_max_base = np.floor((end - start) / min_width)

            if odd:
                n_min = n_min_base + n_min_base % 2
                n_max = n_max_base - n_max_base % 2
            else:
                n_min = n_min_base + 1 - n_min_base % 2
                n_max = n_max_base - 1 + n_max_base % 2

            logger.debug(
                f"odd {odd} n_min {n_min_base} {n_min}" f" n_max {n_max_base} {n_max}"
            )

            if n_min > n_max:
                # Can't put even one, continue side 0
                sb = sb.choose_interside(sides[0].id, geom1d)
                logger.debug(f"None possible choose side 0: f{sides[0].id}")
                continue

            preferred_min = 5

            if n_min < preferred_min:
                n_min = min(preferred_min, n_max)

            n = rng.randint(n_min, n_max + 1)

            locations = random_min_max(n, min_width, max_width, end - start, rng)

            cur_side = 1 - on_side[0]

            side_fingers = [Geom1D.empty(), Geom1D.empty()]

            for span_0, span_1 in locations:
                span = Geom1D([(start + F(span_0), start + F(span_1))])
                logger.debug("Choose %s %s", sides[cur_side].id, span)
                side_fingers[cur_side] |= span
                cur_side = 1 - cur_side

            for side in range(2):
                sb = sb.choose_interside(sides[side].id, side_fingers[side])
                # sb.check()
            # sb.check()
    return sb


def choose_all_ok(sb, params=None):
    """Choose all "ok" areas to be in.

    To be run after all conflicts resolved.
    """
    sp = sb.sheetplex
    for sheet in sp.sheets.values():
        sb = sb.choose_sheet_area(sheet.id, sb.sheet_ok[sheet.id])
    return sb


def clean_up_thin_ok_parts(sb, params=None):
    sp = sb.sheetplex
    eps_pos = F(1, 100)  # noqa
    eps_neg = F(1, 2)
    sheet_ok = sb.sheet_ok
    for sheet in sp.sheets.values():
        logger.debug(sheet.id)
        chosen = sb.sheet_chosen[sheet.id]
        ok = sb.sheet_ok[sheet.id]
        ok = ok.buffer(-eps_neg).buffer(eps_neg)
        # ok = ok.buffer(eps_pos).buffer(-eps_pos - eps_neg).buffer(eps_neg)
        # Shapely blows up without the buffer
        ok = ok | chosen  # .buffer(-F(1, 1000))
        sheet_ok = sheet_ok.set(sheet.id, ok)
    return sb.set("sheet_ok", sheet_ok)


def select_all(sb, params=None):
    assert 0 == 1, "Doesn't work right for 4-crossings!"
    return sb.set("sheet_chosen", sb["sheet_ok"])
