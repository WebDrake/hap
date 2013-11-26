A class-based attempt at a 2nd-generation random module for D.

The choice of classes rather than structs reflects the fundamental
requirement that random number generators operate as reference types:
the lack of this leads to the current unfortunate situation in
std.random where objects like RandomCover or RandomSample that wrap
an underlying RNG can only copy it by value, which may in turn lead
to many unintended correlations in supposedly random output.

The broad plan is for std.random2 to be a package comprising the
following modules, which are selected by analogy to the corresponding
functionality in the C++11 standard:

   std.random2.generator -- uniform (pseudo-)random number generators

   std.random2.distribution -- random distributions (e.g. uniform,
                               normal, pareto, etc.)

   std.random2.device -- "true" sources of randomness

   std.random2.adaptor -- objects that wrap RNGs and transform them
                          into other forms of randomness
