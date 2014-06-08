hap.random
==========

A D library for random number generation and related functionality such
as sampling, shuffling and so on.

`hap.random` is intended as a replacement for Phobos' `std.random`, and
addresses a number of issues encountered in that module.  In particular,
`hap.random` implements random number generators and related entities
as reference types rather than value types.  It is however largely
derivative of `std.random` and to the greatest extent possible
implements the same API, with some functional additions, notably the
random distribution ranges.

The library name, `hap`, is Welsh for "random" or "chance". :-)


Using
-----

To import all standard functionality from the library, use

```
import hap.random;
```

The individual package modules `hap.random.generator`, `hap.random.distribution`,
`hap.random.adaptor` and `hap.random.traits` can alternatively be imported
individually depending on need.

The experimental module `hap.random.device` must be imported individually
and should be used with caution; its API may be subject to change.

Migration instructions for `std.random` users are provided in the documentation
at http://code.braingam.es/hap/random/

`hap.random` should be compatible with compilers using the D frontend 2.065 or
later.


Building/installing
-------------------

`hap.random` is a source library, and can be used as a `dub` package
dependency: the dub package is simply called `hap`.

```
dub fetch hap
```

should get you a copy from the `code.dlang.org` repository.  You can add `hap`
as a dependency to your projects by adding the following JSON to the `dub.json`
file:

```JSON
"dependencies": {
    "hap": ">=1.0.0-rc.1"
}
```

In addition, the included Makefile offers some utilities for developers
working directly on the library itself:

```
make unittest
```

will build and run two different unittest executables, `unit` and
`unit-xper`, covering respectively the standard and experimental parts
of the library.

```
make benchmark
```

will build two benchmarking executables: `benchmarknew` benchmarks
`hap.random` itself while `benchmarkold` offers comparable speed tests
of `std.random`.

```
make doc
```

can be used to build documentation.


Contributing
------------

Contributions of code, issue reports or simply suggestions are all welcome.
`hap.random` is released under the Boost License 1.0, so code submissions
must be similarly licensed.
