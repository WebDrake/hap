// Written in the D programming language.
/**
 * Provides access to non-deterministic sources of randomness.  These
 * are implemented as Input Ranges and in many cases are architecture-
 * or OS-dependent.
 *
 * As with the pseudo-random number generators provided elsewhere in
 * this package, all random devices are implemented as final classes
 * to ensure reference semantics.
 *
 * Copyright: © 2014 Joseph Rushton Wakeling
 *
 * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0).
 *
 * Authors: $(WEB braingam.es, Joseph Rushton Wakeling)
 *
 * Source: $(PHOBOSSRC std/random2/_device.d)
 */
module std.random2.device;

import std.random2.traits;

import std.range, std.traits;

unittest
{
    import std.stdio;
    writeln("Imported std.random2.device!");
}
