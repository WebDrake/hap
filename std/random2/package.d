// Written in the D programming language.

/**
 * Implements random number generation and related functionality such as
 * sampling, shuffling and so on.  Some of the functions and entities
 * provided here are semantically equivalent or similar to those in the
 * $(D $(LESS)_random$(GREATER)) header in the C++ Standard Template
 * Library.
 *
 * The functionality implemented by this package is divided into four
 * different groups:
 *
 * $(UL
 *   $(LI $(B random number generators), deterministic pseudo-random algorithms;)
 *   $(LI $(B random distributions);)
 *   $(LI $(B random adaptors), such as shuffling or sampling objects;)
 *   $(LI $(B traits) related to random number generation.)
 * )
 *
 * The std._random2 package will import all of the above functionality.
 * Alternatively, individual modules can be imported as required.
 *
 * Experimental functionality, not fully developed but available as a
 * technology preview, includes:
 *
 * $(UL
 *   $(LI $(B random devices), non-deterministic sources of randomness;)
 * )
 *
 * This functionality will not be imported as part of the main std.random2
 * package but can be imported via the individual modules.
 *
 * This package is intended as a replacement for the existing Phobos
 * std.random, and corrects a number of issues found in that module.  In
 * particular, it implements random number generators and related entities as
 * reference types rather than value types.  It is however largely derivative
 * of std.random and to the greatest extent possible implements the same API,
 * with some functional additions, notably the random distribution range types.
 *
 * Bear in mind that non-reference-type RNGs used in conjunction with this
 * package will almost certainly generate erroneous results.  In particular
 * this package should not be used together with std.random itself.
 *
 * Copyright: © 2008-2011 Andrei Alexandrescu,
 *              2013      Nils Boßung,
 *              2013      Chris Cain,
 *              2011      Masahiro Nakagawa,
 *              2012-2014 Joseph Rushton Wakeling
 *
 * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0).
 *
 * Authors: $(WEB erdani.org, Andrei Alexandrescu),
 *          Nils Boßung,
 *          Chris Cain,
 *          Masahiro Nakagawa,
 *          $(WEB braingam.es, Joseph Rushton Wakeling)
 *
 * Source: $(PHOBOSSRC std/_random2/package.d)
 */
module std.random2;

public import std.random2.adaptor;
public import std.random2.device;
public import std.random2.distribution;
public import std.random2.generator;
public import std.random2.traits;
