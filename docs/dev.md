# Internals for developers

## Coding standard

* PEP8 (autopep8).

* flake8 must pass with no warnings.

* Prefer immutable objects wherever appropriate.

* Use exact rational numbers wherever appropriate.
  Do approximate operations only at points where it does not
  hurt the exactness of the results.
  The reason for this is that using exact 2D geometry
  gets rid of a lot of possible problem cases due to
  e.g. triangle flipping due to rounding errors.

## Plycutter process

There are three different stages to plycutter:

* Reading the .STL file into a `SheetPlex`.

* Creating a "free" `SheetBuild` 
  object and incrementally making heuristic choices
  within it.

* Post-processing and writing out the .dxf files.

### SheetPlex

A `SheetPlex` is a mostly policy-free representation
of the input model, after finding where the possible
sheets in the input model are.

It contains slices through the input model to provide
information for the heuristics as well as the relationships
between the sheets.

Intersections between sheets are represented,
as well as the projection of the intersection to either of
the sheets involved (this is important for the heuristics)

### SheetBuild

A `SheetBuild` is the "blackboard" object used by the heuristic
routines. It starts life as a very open description
of the situation
("there could be material here on this sheet and here")
and as the heuristics progress, they make decisions
and convert some possible points to certainty
and some to impossibility, thereby creating the joint
between two sheets.

A SheetBuild is a persistent map implemented
using `pyrsistent`. This way, none of the functions *modify*
anything but only return a new version.

This makes debugging the heuristics a *lot* easier since
it is easy to save the state at each point in the heuristics
chain and rerun a particular step with changed code,
with confidence that things are as they should be.

### Heuristics

The heuristics start by looking at the proposed joints
to see which parts are clearly not meant to be implemented
by a particular sheet (just slicing the 3D model produces
surprising things here that the heuristics mostly remove).

After this, the heuristics look at multi-intersections, i.e.,
the intersections of more than two sheets since those regions
require special care.

Currently, one sheet is chosen for all of such an area
(this is a feature where improvements are still needed).

After this, the heuristics look at the two-sheet intersections
and generate fingers there.

The heuristics are in the package `plycutter.heuristics`
and their driver is in `plycutter.canned`.

### Writeout

Currently, the only postprocessing is the kerf compensation
by a fixed amount. More could be added here, such as
dog-bone compensation (though that might also belong
in an earlier stage, this needs to be figured out when adding it)

## Libraries

Plycutter includes some libraries to provide a good API
for implementing the heuristics.

### Geometry

The main library components are the `Geom1D` and
`Geom2D` classes.
They implement exact polygonal sets in 1D and 2D space,
in a way that disallows zero-measure object droppings
or cavities.
This is the same as 'regularized boolean set-operations' in
CGAL.

#### Geom1D

#### Geom2D

The implementation of `Geom2D`
is the largest part of plycutter currently; replacing
this with a better, external implementation would be wonderful.



