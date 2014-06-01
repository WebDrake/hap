// Written in the D programming language.

/**
 * Implements compile-time checks for different features of random
 * number generating code.
 *
 * Copyright: © 2008-2011 Andrei Alexandrescu,
 *              2013-2014 Joseph Rushton Wakeling
 *
 * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0).
 *
 * Authors: $(WEB erdani.org, Andrei Alexandrescu),
 *          $(WEB braingam.es, Joseph Rushton Wakeling)
 *
 * Source: $(PHOBOSSRC hap/random/_traits.d)
 */
module hap.random.traits;

import std.range, std.traits;

/**
 * Test if $(D Range) is a uniform random number generator. The overload
 * taking an $(D ElementType) also makes sure that the RNG generates
 * values of that type.
 *
 * A uniform random number generator has at least the following features:
 * $(UL
 *   $(LI it is an InputRange)
 *   $(LI it has a $(D bool isUniformRandom) field readable in CTFE)
 *   $(LI its element type is a uniform integral type)
 *   $(LI it has $(D min) and $(D max) fields whose type is the same
 *        as the element type)
 * )
 *
 * This quite strict definition follows that in the C++11 standard.
 * Note that the template is unable to enforce some required features,
 * such as the requirement that the RNG's values must be drawn from the
 * $(I closed) interval $(D [min, max]).
 */
template isUniformRNG(Range, ElementType)
{
    enum bool isUniformRNG = isInputRange!Range &&
        is(typeof(Range.front) == ElementType) &&
        is(typeof(Range.min) == ElementType) &&
        is(typeof(Range.max) == ElementType) &&
        isIntegral!ElementType &&
        isUnsigned!ElementType &&
        is(typeof(
        {
            static assert(Range.isUniformRandom); //tag
        }));
}

/// ditto
template isUniformRNG(Range)
{
    enum bool isUniformRNG =
        is(typeof(
        {
            static assert(isUniformRNG!(Range, typeof(Range.front)));
        }));
}

/**
 * Test if $(D UniformRNG) is a seedable uniform random number generator.
 * The overload taking a $(D SeedType) also makes sure that the generator
 * can be seeded with $(D SeedType).
 *
 * A seedable random-number generator has the following additional features:
 * $(UL
 *   $(LI it has a $(D seed(ElementType)) function)
 * )
 */
template isSeedable(UniformRNG, SeedType)
{
    enum bool isSeedable = isUniformRNG!UniformRNG &&
        is(typeof(
        {
            UniformRNG r = void;    // can define a Rng object
            r.seed(SeedType.init);  // can seed a Rng
        }));
}

///ditto
template isSeedable(UniformRNG)
{
    enum bool isSeedable = isUniformRNG!UniformRNG &&
        is(typeof(
        {
            UniformRNG r = void;           // can define a Rng object
            r.seed(typeof(r.front).init);  // can seed a Rng
        }));
}

unittest
{
    /* Not an RNG because it lacks isUniformRandom,
     * min and max
     */
    struct NoRng
    {
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}
    }
    static assert(!isUniformRNG!(NoRng, uint));
    static assert(!isUniformRNG!(NoRng));
    static assert(!isSeedable!(NoRng, uint));
    static assert(!isSeedable!(NoRng));
    NoRng noRng;
    noRng.popFront();
    assert(!noRng.empty);
    assert(!noRng.front);

    /* Not an RNG because isUniformRandom is false,
     * and it lacks min and max
     */
    struct NoRng2
    {
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}

        enum isUniformRandom = false;
    }
    static assert(!isUniformRNG!(NoRng2, uint));
    static assert(!isUniformRNG!(NoRng2));
    static assert(!isSeedable!(NoRng2, uint));
    static assert(!isSeedable!(NoRng2));
    NoRng2 noRng2;
    noRng2.popFront();
    assert(!noRng2.empty);
    assert(!noRng2.front);

    /* Not an RNG because it lacks front, min
     * and max
     */
    struct NoRng3
    {
        @property bool empty() {return false;}
        void popFront() {}

        enum isUniformRandom = true;
    }
    static assert(!isUniformRNG!(NoRng3, uint));
    static assert(!isUniformRNG!(NoRng3));
    static assert(!isSeedable!(NoRng3, uint));
    static assert(!isSeedable!(NoRng3));
    NoRng3 noRng3;
    noRng3.popFront();
    assert(!noRng3.empty);

    // Not an RNG because it lacks min and max
    struct NoRng4
    {
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}

        enum isUniformRandom = true;
    }
    static assert(!isUniformRNG!(NoRng4, uint));
    static assert(!isUniformRNG!(NoRng4));
    static assert(!isSeedable!(NoRng4, uint));
    static assert(!isSeedable!(NoRng4));
    NoRng4 noRng4;
    noRng4.popFront();
    assert(!noRng4.empty);
    assert(!noRng4.front);

    /* Not an RNG because max is different type
     * to front and min
     */
    struct NoRng5
    {
        enum uint min = 0;
        enum ulong max = ulong.max;
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}

        enum isUniformRandom = true;
    }
    static assert(!isUniformRNG!(NoRng5, uint));
    static assert(!isUniformRNG!(NoRng5));
    static assert(!isSeedable!(NoRng5, uint));
    static assert(!isSeedable!(NoRng5));
    NoRng5 noRng5;
    noRng5.popFront();
    assert(!noRng5.empty);
    assert(!noRng5.front);

    /* Not an RNG because its element type is
     * not an unsigned integer
     */
    struct NoRng6
    {
        enum double min = 0;
        enum double max = 23.5;
        @property double front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}

        enum isUniformRandom = true;
    }
    static assert(!isUniformRNG!(NoRng6, double));
    static assert(!isUniformRNG!(NoRng6));
    static assert(!isSeedable!(NoRng6, double));
    static assert(!isSeedable!(NoRng6));
    NoRng6 noRng6;
    noRng6.popFront();
    assert(!noRng6.empty);
    assert(!noRng6.front);

    // Not an RNG because it lacks isUniformRandom
    struct NoRng7
    {
        enum uint min = 0;
        enum uint max = uint.max;
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}
    }
    static assert(!isUniformRNG!(NoRng7, uint));
    static assert(!isUniformRNG!(NoRng7));
    static assert(!isSeedable!(NoRng7, uint));
    static assert(!isSeedable!(NoRng7));
    NoRng7 noRng7;
    noRng7.popFront();
    assert(!noRng7.empty);
    assert(!noRng7.front);

    // Not an RNG because isUniformRandom is false
    struct NoRng8
    {
        enum uint min = 0;
        enum uint max = uint.max;
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}

        enum isUniformRandom = false;
    }
    static assert(!isUniformRNG!(NoRng8, uint));
    static assert(!isUniformRNG!(NoRng8));
    static assert(!isSeedable!(NoRng8, uint));
    static assert(!isSeedable!(NoRng8));
    NoRng8 noRng8;
    noRng8.popFront();
    assert(!noRng8.empty);
    assert(!noRng8.front);

    // Valid RNG, but not seedable
    struct ValidRng
    {
        enum uint min = 0;
        enum uint max = uint.max;
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}

        enum isUniformRandom = true;
    }
    static assert(isUniformRNG!(ValidRng, uint));
    static assert(isUniformRNG!(ValidRng));
    static assert(!isSeedable!(ValidRng, uint));
    static assert(!isSeedable!(ValidRng));
    ValidRng validRng;
    validRng.popFront();
    assert(!validRng.empty);
    assert(!validRng.front);

    // Valid and seedable RNG
    struct SeedRng
    {
        enum ulong min = 0;
        enum ulong max = ulong.max;
        @property ulong front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}
        void seed(ulong val){}
        enum isUniformRandom = true;
    }
    static assert(isUniformRNG!(SeedRng, ulong));
    static assert(isUniformRNG!(SeedRng));
    static assert(isSeedable!(SeedRng, ulong));
    static assert(isSeedable!(SeedRng));
    SeedRng seedRng;
    seedRng.seed(123456789uL);
    seedRng.popFront();
    assert(!seedRng.empty);
    assert(!seedRng.front);
}


/**
 * Test if $(D RandomDist) is a random distribution. The overload
 * taking an $(D ElementType) also makes sure that the distribution
 * generates values of that type.
 *
 * A random distribution has at least the following features:
 * $(UL
 *   $(LI it's an InputRange)
 *   $(LI it has a $(D bool isRandomDistribution) field readable
 *        in CTFE)
 * )
 */
template isRandomDistribution(RandomDist, ElementType)
{
    enum bool isRandomDistribution = isInputRange!RandomDist &&
        is(typeof(RandomDist.front) == ElementType) &&
        is(typeof(
        {
            static assert(RandomDist.isRandomDistribution); //tag
        }));
}

/// ditto
template isRandomDistribution(RandomDist)
{
    enum bool isRandomDistribution = isInputRange!RandomDist &&
        is(typeof(
        {
            static assert(RandomDist.isRandomDistribution); //tag
        }));
}

unittest
{
    struct NoDist
    {
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}
    }
    static assert(!isRandomDistribution!(NoDist, uint));
    static assert(!isRandomDistribution!(NoDist));
    NoDist noDist;
    noDist.popFront();
    assert(!noDist.empty);
    assert(!noDist.front);

    struct NoDist2
    {
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}

        enum isRandomDistribution = false;
    }
    static assert(!isRandomDistribution!(NoDist2, uint));
    static assert(!isRandomDistribution!(NoDist2));
    NoDist2 noDist2;
    noDist2.popFront();
    assert(!noDist2.empty);
    assert(!noDist2.front);

    struct NoDist3
    {
        @property bool empty() {return false;}
        void popFront() {}

        enum isRandomDistribution = true;
    }
    static assert(!isRandomDistribution!(NoDist3, uint));
    static assert(!isRandomDistribution!(NoDist3));
    NoDist3 noDist3;
    noDist3.popFront();
    assert(!noDist3.empty);

    struct ValidDist
    {
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}

        enum isRandomDistribution = true;
    }
    static assert(isRandomDistribution!(ValidDist, uint));
    static assert(isRandomDistribution!(ValidDist));
    ValidDist validDist;
    validDist.popFront();
    assert(!validDist.empty);
    assert(!validDist.front);
}
