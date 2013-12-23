module std.random2.traits;

import std.range;

/**
 * Test if $(D RandomGen) is a random-number generator. The overload
 * taking an $(D ElementType) also makes sure that the Rng generates
 * values of that type.
 *
 * A random-number generator has at least the following features:
 * $(UL
 *   $(LI it's an InputRange)
 *   $(LI it has a $(D bool isUniformRandom) field readable in CTFE)
 * )
 */
template isUniformRNG(RandomGen, ElementType)
{
    enum bool isUniformRNG = isInputRange!RandomGen &&
        is(typeof(RandomGen.front) == ElementType) &&
        is(typeof(
        {
            static assert(RandomGen.isUniformRandom); //tag
        }));
}

/// ditto
template isUniformRNG(RandomGen)
{
    enum bool isUniformRNG = isInputRange!RandomGen &&
        is(typeof(
        {
            static assert(RandomGen.isUniformRandom); //tag
        }));
}

/**
 * Test if $(D RandomGen) is a seedable uniform random number generator.
 * The overload taking a $(D SeedType) also makes sure that the generator
 * can be seeded with $(D SeedType).
 *
 * A seedable random-number generator has the following additional features:
 * $(UL
 *   $(LI it has a $(D seed(ElementType)) function)
 * )
 */
template isSeedable(RandomGen, SeedType)
{
    enum bool isSeedable = isUniformRNG!RandomGen &&
        is(typeof(
        {
            RandomGen r = void;     // can define a Rng object
            r.seed(SeedType.init);  // can seed a Rng
        }));
}

///ditto
template isSeedable(RandomGen)
{
    enum bool isSeedable = isUniformRNG!RandomGen &&
        is(typeof(
        {
            RandomGen r = void;            // can define a Rng object
            r.seed(typeof(r.front).init);  // can seed a Rng
        }));
}

unittest
{
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

    struct ValidRng
    {
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

    struct SeedRng
    {
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}
        void seed(uint val){}
        enum isUniformRandom = true;
    }
    static assert(isUniformRNG!(SeedRng, uint));
    static assert(isUniformRNG!(SeedRng));
    static assert(isSeedable!(SeedRng, uint));
    static assert(isSeedable!(SeedRng));
    SeedRng seedRng;
    seedRng.seed(123456789U);
    seedRng.popFront();
    assert(!seedRng.empty);
    assert(!seedRng.front);
}

