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

import matchpy

import hypothesis as hyp
import hypothesis.stateful
import hypothesis.strategies as hys

from ..geom1d import Geom1D

from ...plytypes import Fraction

# hyp.settings.register_profile('long', max_examples=2000)
# hyp.settings.load_profile('long')


@hys.composite
def ply_fractions(draw):
    return Fraction(draw(hys.fractions()))


@hys.composite
def random_transforms(draw):
    c = [draw(ply_fractions()) for _ in range(2)]
    hyp.assume(c[0] != 0)
    return c


def tr(transform, pt1):
    if transform is None:
        return pt1
    return transform[0] * pt1 + transform[1]


F = Fraction


def test_basics():
    # inf = math.inf
    assert Geom1D([[0, 1], [1, 2]]).intervals == [(0, 2)]
    assert Geom1D([[0, 1], [2, 3]]).intervals == [(0, 1), (2, 3)]

    assert (Geom1D([[0, 1]]) & Geom1D(
        [[F(0.5), F(1.5)]])).intervals == [(F(0.5), 1)]
    assert (Geom1D([[0, 1]]) | Geom1D(
        [[F(0.5), F(1.5)]])).intervals == [(0, F(1.5))]

    assert (Geom1D.empty() & ~Geom1D([[-200, 0]])).is_empty()

# A test version of the 2D tester...


def rule(pat, call):
    return matchpy.ReplacementRule(
        matchpy.Pattern(pat), call)


_a = matchpy.Wildcard.dot('a')
_amount = matchpy.Wildcard.dot('amount')
_b = matchpy.Wildcard.dot('b')
_c = matchpy.Wildcard.dot('c')
op_and = matchpy.Operation.new('op_and', matchpy.Arity.binary)
op_or = matchpy.Operation.new('op_or', matchpy.Arity.binary)
op_sub = matchpy.Operation.new('op_sub', matchpy.Arity.binary)
op_buffer = matchpy.Operation.new('op_buffer', matchpy.Arity.binary)


class SingleInterval(matchpy.Symbol):
    def __init__(self, vertices):
        super().__init__(str(vertices))
        self.vertices = vertices


equivalence_rules = [
    rule(op_and(_a, _b),
         lambda a, b: op_and(b, a)),
    rule(op_and(_a, op_and(_b, _c)),
         lambda a, b, c: op_and(op_and(a, b), c)),

    rule(op_or(_a, _b),
         lambda a, b: op_or(b, a)),
    rule(op_or(_a, op_or(_b, _c)),
         lambda a, b, c: op_or(op_or(a, b), c)),

    rule(op_and(_a, op_or(_b, _c)),
         lambda a, b, c: op_or(op_and(a, b), op_and(a, c))),

    rule(op_sub(_a, op_or(_b, _c)),
         lambda a, b, c: op_and(op_sub(a, b), op_sub(a, c))),
    rule(op_sub(op_or(_a, _b), _c),
         lambda a, b, c: op_or(op_sub(a, c), op_sub(b, c))),
]


def is_equivalent(g0, g0tr, g1, g1tr):
    pts = []
    eps = Fraction(1, 1_000_000_000)
    for intervals in [g0.intervals, g1.intervals]:
        for interval in intervals:
            pts.append(interval[0])
            pts.append(interval[0] + eps)
            pts.append(interval[0] - eps)
            pts.append(interval[1])
            pts.append(interval[1] + eps)
            pts.append(interval[1] - eps)
            pts.append((interval[0] + interval[1]) / 2)
    for pt in pts:
        r0 = g0.locate(tr(g0tr, pt))
        r1 = g1.locate(tr(g1tr, pt))
        if r0 == r1:
            continue
        if r0 >= 0 and r1 >= 0:
            continue
        return False
    return True


class GeomExpressionTree(hyp.stateful.RuleBasedStateMachine):
    """Assert invariants on generated expression trees.

    We generate a bunch of representations (i.e.
    affine transformations in which we run the same ops),
    and a bunch of expression trees which we evaluate
    in all transforms and on which we run various equivalence
    transforms.

    We can then assert that the things that should be equivalent
    really are
    """

    # Names of the expressions used
    expr_names = hyp.stateful.Bundle('expr_names')

    # Indices of the representations
    repr_inds = hyp.stateful.Bundle('repr_inds')

    def evaluate(self, expr):
        # print('Eval', expr)
        cached = self.cache.get(expr, None)
        if cached is not None:
            return cached

        eva = self.evaluate
        eval_rules = [
            rule(op_and(_a, _b),
                 lambda a, b: eva(a) & eva(b)),
            rule(op_or(_a, _b),
                 lambda a, b: eva(a) | eva(b)),
            rule(op_sub(_a, _b),
                 lambda a, b: eva(a) - eva(b)),
            rule(matchpy.Wildcard.symbol('interval', SingleInterval),
                 lambda interval:
                 Geom1D([[interval.vertices[0],
                          interval.vertices[1]]])),
        ]
        for eval_rule in eval_rules:
            matches = list(matchpy.match(expr, eval_rule.pattern))
            if matches:
                assert len(matches) == 1
                result = eval_rule.replacement(**matches[0])
                self.cache[expr] = result
                return result
        assert 0 == 1, expr

    def check(self, name):
        expr = self.exprs[name]
        # print('check expr', expr)
        expr = [self.evaluate(e) for e in expr]
        for i in range(1, len(expr)):
            assert is_equivalent(
                expr[0], self.reprs[0],
                expr[i], self.reprs[i]), i

    @hyp.stateful.initialize(
        target=repr_inds,
        transforms=hys.lists(elements=random_transforms()))
    def init_reprs(self, transforms):
        self.exprs = {}
        self.cache = {}
        self.reprs = [[1, 0]] + transforms
        return hyp.stateful.multiple(
            *list(range(len(self.reprs))))

    @hyp.stateful.rule(target=expr_names,
                       name=hys.text("abcdefghijklmnopqrstuvwxyz0123456789",
                                     ),
                       element=hys.tuples(ply_fractions(), ply_fractions()))
    def add_single(self, name, element):
        res = []
        for rep in self.reprs:
            res.append(
                SingleInterval(tuple(sorted(
                    [tr(rep, coord) for coord in element]))))
        self.exprs[name] = res
        self.check(name)
        return name

    @hyp.stateful.rule(
        dst=expr_names, a=expr_names, b=expr_names,
        op=hys.sampled_from([
            op_and,
            op_or, op_sub
        ]))
    def binary_op(self, dst, op, a, b):
        self.exprs[dst] = [op(*params)
                           for params in zip(self.exprs[a], self.exprs[b])]
        self.check(dst)

    @hyp.stateful.rule(
        name=expr_names,
        idx=repr_inds,
        rules=hys.lists(
            elements=hys.sampled_from(equivalence_rules)))
    def equivalence_transform(self, name, idx, rules):
        expr = [e for e in self.exprs[name]]
        for rule in rules:
            expr[idx] = matchpy.replace_all(expr[idx],
                                            [rule], max_count=1)

        self.exprs[name] = expr
        self.check(name)


TestGeomExpressionTree = GeomExpressionTree.TestCase
