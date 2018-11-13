Dune-Vof
========

`dune-vof` is a Dune module that implements a volume-of-fluid method within the Dune framework.

It enables interface reconstructions for a piecewise constant color function up to three dimensions and the sharp tracking of the interface movement.

You can have a look at some [videos of results](https://av.tib.eu/series/350/volume+of+fluid+examples) computed with this module.

Implemented interface reconstructions are:

| Reconstruction | <abbr title="Supported grid types">Grid</abbr>  | <abbr title="Order of error between reconstructed and exact interface">Order</abbr> | <abbr title="Computational effort">Cost</abbr> |
|---------------------------------|-----------|-------|------|
| Young's Least Squares           | all       | 1     | Low  |
| Iterative Swartz reconstruction | all       | 2     | High |
| Height function technique       | cartesian | 2     | Low  |


There is also a grid interface called `interfacegrid` included which allows to iterate over the reconstructed surface like a Dune grid.

Moreover there are methods to compute the curvature of the reconstructed interface.


Usage
=====

You can install and use this package as every other external dune module.  
As starting point have a look at the `/example/vof.cc`.
