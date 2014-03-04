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
 *   $(LI $(B random number generators);)
 *   $(LI $(B random distributions);)
 *   $(LI $(B random adaptors), such as shuffling or sampling objects;)
 *   $(LI $(B traits) related to random number generation.)
 * )
 *
 * Copyright: © 2008-2011 Andrei Alexandrescu,
 *              2011      Masahiro Nakagawa,
 *              2012-2014 Joseph Rushton Wakeling
 *
 * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0).
 *
 * Authors: $(WEB erdani.org, Andrei Alexandrescu),
 *          Masahiro Nakagawa (Xorshift random generator),
 *          $(WEB braingam.es, Joseph Rushton Wakeling)
 *
 * Source: $(PHOBOSSRC std/_random2/package.d)
 */
module std.random2;

public import std.random2.adaptor;
public import std.random2.distribution;
public import std.random2.generator;
public import std.random2.traits;
