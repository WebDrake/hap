A class-based attempt at a 2nd-generation random module for D.

The choice of classes rather than structs reflects the fundamental
requirement that random number generators operate as reference types:
the lack of this leads to the current unfortunate situation in
std.random where objects like RandomCover or RandomSample that wrap
an underlying RNG can only copy it by value, which may in turn lead
to many unintended correlations in supposedly random output.

The broad plan is for hap.random to be a package comprising the
following modules, most of which are selected by analogy to the
corresponding functionality in the C++11 standard:

   hap.random.traits -- compile-time checks and templates for working
                        with random number generators and related code

   hap.random.generator -- uniform (pseudo-)random number generators

   hap.random.distribution -- random distributions (e.g. uniform,
                              normal, pareto, etc.)

   hap.random.device -- "true" sources of randomness

   hap.random.adaptor -- objects that wrap RNGs and transform them
                         into other forms of randomness
