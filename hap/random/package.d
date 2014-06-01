// Written in the D programming language.

/**
 * Implements _random number generation and related functionality such as
 * sampling, shuffling and so on.  Some of the functions and entities
 * provided here are semantically equivalent or similar to those in the
 * $(D $(LESS)_random$(GREATER)) header in the C++ Standard Template
 * Library.
 *
 * $(D hap._random) is intended as a replacement for Phobos' $(D std._random),
 * and addresses a number of issues encountered in that module.  In particular,
 * $(D hap._random) implements _random number generators and related entities as
 * reference types rather than value types.  It is however largely derivative
 * of $(D std._random) and to the greatest extent possible implements the same
 * API, with some functional additions, notably the _random distribution ranges.
 *
 * The library name, $(I hap), is Welsh for "_random" or "chance". :-)
 *
 * The functionality implemented by this package is divided into four
 * different groups:
 *
 * $(UL
 *   $(LI $(B _random number generators), deterministic pseudo-random algorithms;)
 *   $(LI $(B _random distributions);)
 *   $(LI $(B _random adaptors), such as shuffling or sampling objects;)
 *   $(LI $(B traits) related to _random number generation.)
 * )
 *
 * The $(D hap._random) package will import all of the above functionality.
 * Alternatively, individual modules can be imported as required.
 *
 * Experimental functionality, not fully developed but available as a
 * technology preview, includes:
 *
 * $(UL
 *   $(LI $(B _random devices), non-deterministic sources of randomness.)
 * )
 *
 * This functionality will not be imported as part of the main $(D hap._random)
 * package but can be imported via the individual modules.  It should be used
 * with some caution, as its API may change significantly in future releases.
 *
 * Migration: To use $(D hap._random) instead of $(D std._random) it should
 * suffice to:
 *
 * $(UL
 *   $(LI replace $(D import std._random) with $(D import hap._random))
 *   $(LI insert $(D new) before every instantiation of an RNG.  For example,
 *     instead of "$(D auto rng = Random(unpredictableSeed))" put instead,
 *     "$(D auto rng = new Random(unpredictableSeed))".)
 *   $(LI optionally, the functions $(D randomCover), $(D randomSample) and
 *     $(D randomShuffle) may be replaced by $(D cover), $(D sample) and
 *     $(D shuffle); however, this is optional, not required.)
 * )
 *
 * Note: Currently only one function is known to produce different behaviour to
 * its counterpart in $(D std._random): $(D hap._random.distribution.dice) uses
 * a different algorithm to $(D std._random.dice).
 *
 * Warning: Bear in mind that non-reference-type RNGs used in conjunction with
 * this package will almost certainly generate erroneous results.  In particular
 * this package should not be used together with $(D std._random) itself.
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
 * Source: $(HAPSRC hap/_random/package.d)
 */
module hap.random;

public import hap.random.adaptor;
public import hap.random.distribution;
public import hap.random.generator;
public import hap.random.traits;
