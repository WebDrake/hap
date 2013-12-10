module std.random2.generator;

import std.algorithm, std.conv, std.exception, std.range, std.traits, std.typetuple, core.thread;
import std.string : format;

unittest
{
    import std.stdio;
    writeln("std.random2.generator has been imported.");
}

/// Tuple of all uniform RNGs defined in this module.
alias UniformRNGTypes = TypeTuple!(MinstdRand0, MinstdRand, Mt11213b, Mt19937, Mt19937_64);

/// Default RNG type recommended for general use.
alias Random = Mt19937;

ref Random rndGen() @property
{
    struct RndGen
    {
        Random gen = new Random();
    }

    static RndGen result;
    static bool initialized;

    if (!initialized)
    {
        static if (isSeedable!(Random, typeof(repeat(0).map!((a) => unpredictableSeed))))
        {
            result.gen.seed(repeat(0).map!((a) => unpredictableSeed));
        }
        else
        {
            result.gen.seed(unpredictableSeed);
        }

        initialized = true;
    }

    return result.gen;
}

unittest
{
    assert(isUniformRNG!(typeof(rndGen)));
    auto a = rndGen;
    assert(a.front == rndGen.front);
}

// General unittests that all uniform RNGs should pass
unittest
{
    foreach (RandomGen; UniformRNGTypes)
    {
        assert(isUniformRNG!RandomGen);
        assert(isSeedable!RandomGen);

        import std.stdio;
        writeln(RandomGen.stringof);
        writeln("\tType: ", typeof(RandomGen.front).stringof);
        writeln("\tmin = ", RandomGen.min);
        writeln("\tmax = ", RandomGen.max);

        // Ensure that popFront() actually changes the RNG state
        typeof(RandomGen.front) a, b;
        {
            auto gen = new RandomGen;
            a = gen.front;
        }
        {
            auto gen = new RandomGen;
            gen.popFront();
            b = gen.front;
        }
        assert(a != b);

        static if (isForwardRange!RandomGen)
        {
            auto gen1 = new RandomGen(unpredictableSeed);
            auto gen2 = gen1.save;
            assert(gen1 == gen2);
            assert(gen1 !is gen2);
            assert(gen1.take(100).array() == gen2.take(100).array());
            auto gen3 = gen1;
            gen1.popFrontN(9999);
            assert(gen3 == gen1);
            assert(gen3 is gen1);
            assert(gen3.front == gen1.front);
            assert(gen2 != gen1);
            assert(gen2 !is gen1);
            assert(gen2.front != gen1.front);
        }
    }
}

/**
 * Test if Rng is a random-number generator. The overload
 * taking a ElementType also makes sure that the Rng generates
 * values of that type.
 *
 * A random-number generator has at least the following features:
 * $(UL
 *   $(LI it's an InputRange)
 *   $(LI it has a 'bool isUniformRandom' field readable in CTFE)
 * )
 */
template isUniformRNG(Rng, ElementType)
{
    enum bool isUniformRNG = isInputRange!Rng &&
        is(typeof(Rng.front) == ElementType) &&
        is(typeof(
        {
            static assert(Rng.isUniformRandom); //tag
        }));
}

/**
 * ditto
 */
template isUniformRNG(Rng)
{
    enum bool isUniformRNG = isInputRange!Rng &&
        is(typeof(
        {
            static assert(Rng.isUniformRandom); //tag
        }));
}

/**
 * Test if Rng is seedable. The overload
 * taking a SeedType also makes sure that the Rng can be seeded with SeedType.
 *
 * A seedable random-number generator has the following additional features:
 * $(UL
 *   $(LI it has a 'seed(ElementType)' function)
 * )
 */
template isSeedable(Rng, SeedType)
{
    enum bool isSeedable = isUniformRNG!(Rng) &&
        is(typeof(
        {
            Rng r = void;              // can define a Rng object
            r.seed(SeedType.init);     // can seed a Rng
        }));
}

///ditto
template isSeedable(Rng)
{
    enum bool isSeedable = isUniformRNG!Rng &&
        is(typeof(
        {
            Rng r = void;                     // can define a Rng object
            r.seed(typeof(r.front).init);     // can seed a Rng
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

    struct validRng
    {
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}

        enum isUniformRandom = true;
    }
    static assert(isUniformRNG!(validRng, uint));
    static assert(isUniformRNG!(validRng));
    static assert(!isSeedable!(validRng, uint));
    static assert(!isSeedable!(validRng));

    struct seedRng
    {
        @property uint front() {return 0;}
        @property bool empty() {return false;}
        void popFront() {}
        void seed(uint val){}
        enum isUniformRandom = true;
    }
    static assert(isUniformRNG!(seedRng, uint));
    static assert(isUniformRNG!(seedRng));
    static assert(isSeedable!(seedRng, uint));
    static assert(isSeedable!(seedRng));
}

final class LinearCongruentialEngine(UIntType,
                                     UIntType a, UIntType c, UIntType m)
    if (isUnsigned!UIntType && isIntegral!UIntType)
{
  private:
    UIntType _x = m ? (a + c) % m : (a + c);

    static ulong primeFactorsOnly(ulong n) @safe nothrow pure
    {
        ulong result = 1;
        ulong iter = 2;
        for(; n >= iter * iter; iter += 2 - (iter == 2))
        {
            if (n % iter) continue;
            result *= iter;
            do
            {
                n /= iter;
            }
            while (n % iter == 0);
        }
        return result * n;
    }

    unittest
    {
        static assert(primeFactorsOnly(100) == 10);
        static assert(primeFactorsOnly(11) == 11);
        static assert(primeFactorsOnly(7 * 7 * 7 * 11 * 15 * 11) == 3 * 5 * 7 * 11);
        static assert(primeFactorsOnly(129 * 2) == 129 * 2);
    }

    static bool properParameters(ulong m, ulong a, ulong c) @safe nothrow pure
    {
        if (m == 0)
        {
            static if (is(UIntType == uint))
            {
                // assume m is uint.max + 1
                m = 1uL << 32;
            }
            else
            {
                return false;
            }
        }

        import std.numeric : gcd;

        // Bounds checking
        if (a == 0 || a >= m || c >= m) return false;
        // c and m are relatively prime
        if (c > 0 && gcd(c, m) != 1) return false;
        // a - 1 is divisible by all prime factors of m
        if ((a - 1) % primeFactorsOnly(m)) return false;
        // if a - 1 is a multiple of 4, then m is a multiple of 4 too
        if ((a - 1) % 4 == 0 && m % 4) return false;

        // Passed all tests
        return true;
    }

    static assert(c == 0 || properParameters(m, a, c), format("Incorrect instantiation of ", typeof(this).stringof));

  public:
    enum bool isUniformRandom = true;
    enum UIntType min = (c == 0 ? 1 : 0);
    enum UIntType max = m - 1;

    enum UIntType multiplier = a;
    enum UIntType increment = c;
    enum UIntType modulus = m;

    static assert(m == 0 || a < m);
    static assert(m == 0 || c < m);
    static assert(m == 0 || (cast(ulong) a * (m - 1) + c) % m == (c < a ? c - a + m : c - a));

    this()
    {
        // Default seed already set to m ? (a + c) % m : (a + c)
    }

    this(in UIntType x0)
    {
        seed(x0);
    }

    void seed(in UIntType x0 = 1) @safe pure
    {
        static if (c == 0)
        {
            enforce(x0, format("Invalid (zero) seed for %s", typeof(this).stringof));
        }
        _x = m ? (x0 % m) : x0;
        popFront();
    }

    // ----- Range primitives -------------------------------------------------

    /// Always $(D false) (random generators are infinite ranges).
    enum bool empty = false;

    /// Returns the current pseudo-random value.
    UIntType front() @property @safe const nothrow pure
    {
        return _x;
    }

    void popFront() @safe nothrow pure
    {
        static if (m)
        {
            static if (is(UInttype == uint) && m == uint.max)
            {
                immutable ulong x = (cast(ulong) a * _x + c),
                                v = x >> 32,
                                w = x & uint.max;
                immutable y = cast(uint)(v + w);
                _x = (y < v || y = uint.max) ? (y + 1) : y;
            }
            else static if (is(UIntType == uint) && m == int.max)
            {
                immutable ulong x = (cast(ulong) a * _x + c),
                                v = x >> 31,
                                w = x & int.max;
                immutable uint y = cast(uint)(v + w);
                _x = (y >= int.max) ? (y - int.max) : y;
            }
            else
            {
                _x = cast(UIntType) ((cast(ulong) a * _x + c) % m);
            }
        }
        else
        {
            _x = a * _x + c;
        }
    }

    typeof(this) save() @property
    {
        auto ret = new typeof(this);
        ret._x = this._x;
        return ret;
    }

    override bool opEquals(Object rhs) @safe const nothrow pure
    {
        auto that = cast(typeof(this)) rhs;

        if (that is null)
        {
            return false;
        }
        else if (this._x != that._x)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}

alias MinstdRand0 = LinearCongruentialEngine!(uint, 16807, 0, 2147483647);
alias MinstdRand = LinearCongruentialEngine!(uint, 48271, 0, 2147483647);

unittest
{
    foreach (LCGen; TypeTuple!(MinstdRand0, MinstdRand))
    {
        assert(isUniformRNG!(LCGen, uint));
        assert(isSeedable!(LCGen, uint));
        assert(isForwardRange!LCGen);
    }

    // The correct numbers are taken from The Database of Integer Sequences
    // http://www.research.att.com/~njas/sequences/eisBTfry00128.txt
    auto checking0 = [
        16807UL,282475249,1622650073,984943658,1144108930,470211272,
        101027544,1457850878,1458777923,2007237709,823564440,1115438165,
        1784484492,74243042,114807987,1137522503,1441282327,16531729,
        823378840,143542612 ];
    auto gen0 = new MinstdRand0;

    foreach (e; checking0)
    {
        assert(gen0.front == e);
        gen0.popFront();
    }
    // Test the 10000th invocation
    // Correct value taken from:
    // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2007/n2461.pdf
    gen0.seed();
    popFrontN(gen0, 9999);
    assert(gen0.front == 1043618065);

    // Test MinstdRand
    auto checking = [48271UL,182605794,1291394886,1914720637,2078669041, 407355683];
    auto gen = new MinstdRand;

    foreach (e; checking)
    {
        assert(gen.front == e);
        gen.popFront();
    }

    // Test the 10000th invocation
    // Correct value taken from:
    // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2007/n2461.pdf
    gen.seed();
    popFrontN(gen, 9999);
    assert(gen.front == 399268537);

    // Check we don't get equality if 2 LC generators of different type are compared.
    gen._x = 1;
    gen0._x = 1;
    assert(gen._x == gen0._x);
    assert(gen !is gen0);
    assert(gen != gen0);
}

final class MersenneTwisterEngine(UIntType,
                                  size_t w, size_t n, size_t m, size_t r,
                                  UIntType a, size_t u, UIntType d, size_t s,
                                  UIntType b, size_t t,
                                  UIntType c, size_t l, UIntType f)
    if (isUnsigned!UIntType)
{
  private:
    UIntType[n] mt;
    UIntType _y;
    size_t mti = size_t.max;   // means mt is not initialized

  public:
    /// Mark this as a uniform RNG
    enum bool isUniformRandom = true;

    /// Parameter for the generator
    enum size_t wordSize = w;
    enum size_t stateSize = n;
    enum size_t shiftSize = m;
    enum size_t maskBits = r;
    enum UIntType xorMask = a;
    enum size_t temperingU = u;
    enum UIntType temperingD = d;
    enum size_t temperingS = s;
    enum UIntType temperingB = b;
    enum size_t temperingT = t;
    enum UIntType temperingC = c;
    enum size_t temperingL = l;
    enum UIntType initializationMultiplier = f;

    /// Smallest generated value (0)
    enum UIntType min = 0;

    /// Largest generated value
    enum UIntType max =
        w == UIntType.sizeof * 8 ? UIntType.max : (to!UIntType(1) << w) - 1;

    /// The default seed value
    enum UIntType defaultSeed = 5489U;

    this()
    {
        seed(this.defaultSeed);
    }

    this(in UIntType value)
    {
        seed(value);
    }

    void seed()(in UIntType value) @safe nothrow pure
    {
        enum UIntType mask = this.max;
        mt[0] = value & mask;
        for (mti = 1; mti < n; ++mti)
        {
            mt[mti] = (f * (mt[mti - 1] ^ (mt[mti - 1] >> (w - 2))) + mti) & mask;
        }
        popFront();
    }

    void seed(Range)(Range range)
    if (isInputRange!Range && is(Unqual!(ElementType!Range) : UIntType))
    {
        size_t j;
        for (j = 0; j < n && !range.empty; ++j, range.popFront())
        {
            mt[j] = range.front;
        }

        mti = n;
        if (range.empty && j < n)
        {
            throw new Exception(format("%s.seed: Input range only provided %s elements, "
                                       "need at least %s.", typeof(this).stringof, j, n));
        }

        popFront();
    }

    // ----- Range primitives -------------------------------------------------

    /// Always $(D false) (random generators are infinite ranges).
    enum bool empty = false;

    /// Returns the current pseudo-random value.
    UIntType front() @property @safe const nothrow pure
    in
    {
        assert(mti < size_t.max);
    }
    body
    {
        return _y;
    }

    /// Advances the generator.
    void popFront() @safe nothrow pure
    in
    {
        assert(mti < size_t.max);
    }
    body
    {
        enum UIntType upperMask = (~(cast(UIntType) 0)) << r;
        enum UIntType lowerMask = ~upperMask;

        enum size_t unrollFactor = 6;
        enum size_t unrollExtra1 = (n - m) % unrollFactor;
        enum size_t unrollExtra2 = (m - 1) % unrollFactor;

        UIntType y = void;

        if (mti >= n)
        {
            foreach (j; 0 .. n - m - unrollExtra1)
            {
                y = (mt[j] & upperMask) | (mt[j + 1] & lowerMask);
                mt[j] = mt[j + m] ^ (y >> 1) ^ ((mt[j + 1] & 1) * a);
            }

            foreach (j; n - m - unrollExtra1 .. n - m)
            {
                y = (mt[j] & upperMask) | (mt[j + 1] & lowerMask);
                mt[j] = mt[j + m] ^ (y >> 1) ^ ((mt[j + 1] & 1) * a);
            }

            foreach (j; n - m .. n - 1 - unrollExtra2)
            {
                y = (mt[j] & upperMask) | (mt[j + 1] & lowerMask);
                mt[j] = mt[j - (n - m)] ^ (y >> 1) ^ ((mt[j + 1] & 1) * a);
            }

            foreach (j; n - 1 - unrollExtra2 .. n - 1)
            {
                y = (mt[j] & upperMask) | (mt[j + 1] & lowerMask);
                mt[j] = mt[j - (n - m)] ^ (y >> 1) ^ ((mt[j + 1] & 1) * a);
            }

            y = (mt[n - 1] & upperMask) | (mt[0] & lowerMask);
            mt[n - 1] = mt[m - 1] ^ (y >> 1) ^ ((mt[0] & 1) * a);

            mti = 0;
        }

        y = mt[mti];
        mti++;
        y ^= ((y >> u) & d);
        y ^= ((y << s) & b);
        y ^= ((y << t) & c);
        y ^= (y >> l);

        _y = y;
    }

    typeof(this) save() @property
    {
        auto ret = new typeof(this);
        ret.mt[] = this.mt[];
        ret._y = this._y;
        ret.mti = this.mti;
        return ret;
    }

    override bool opEquals(Object rhs) @safe const nothrow pure
    {
        auto that = cast(typeof(this)) rhs;
        if (that is null)
        {
            return false;
        }
        else if (this.mt != that.mt || this._y != that._y || this.mti != that.mti)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}

alias Mt11213b =
    MersenneTwisterEngine!(uint,
                           32, 351, 175, 19,
                           0xccab8ee7U, 11, 0xffffffffU, 7,
                           0x31b6ab00U, 15,
                           0xffe50000U, 17, 1812433253U);

alias Mt19937 =
    MersenneTwisterEngine!(uint,
                           32, 624, 397, 31,
                           0x9908b0dfU, 11, 0xffffffffU, 7,
                           0x9d2c5680U, 15,
                           0xefc60000U, 18, 1812433253U);

alias Mt19937_64 =
    MersenneTwisterEngine!(ulong,
                           64, 312, 156, 31,
                           0xb5026f5aa96619e9UL, 29, 0x5555555555555555UL, 17,
                           0x71d67fffeda60000UL, 37,
                           0xfff7eee000000000UL, 43, 6364136223846793005UL);

unittest
{
    auto gen11213b = new Mt11213b;
    gen11213b.popFrontN(9999);
    assert(gen11213b.front == 3809585648U);

    auto gen19937 = new Mt19937;
    gen19937.popFrontN(9999);
    assert(gen19937.front == 4123659995U);

    auto gen19937_64 = new Mt19937_64;
    gen19937_64.popFrontN(9999);
    assert(gen19937_64.front == 9981545732273789042UL);
}

unittest
{
    foreach (MtGen; TypeTuple!(Mt11213b, Mt19937, Mt19937_64))
    {
        assert(isSeedable!(MtGen, typeof(MtGen.front)));
        assert(isForwardRange!MtGen);

        auto gen = new MtGen;

        // range too small to seed the generator
        assertThrown(gen.seed(repeat(0, gen.stateSize - 1).map!((a) => (a + 1))));

        // range the absolute minimum size necessary
        gen.seed(repeat(0, gen.stateSize).map!((a) => (a + 1)));

        // infinite range
        gen.seed(repeat(0).map!((a) => (a + 1)));
    }
}

uint unpredictableSeed() @property
{
    static MinstdRand0 rand = null;

    if (rand is null)
    {
        rand = new MinstdRand0;
        uint threadID = cast(uint) cast(void*) Thread.getThis();
        rand.seed((getpid() + threadID) ^ cast(uint) TickDuration.currSystemTick.length);
    }
    rand.popFront();
    return cast(uint) TickDuration.currSystemTick.length ^ rand.front;
}

unittest
{
    auto a = unpredictableSeed;
    assert(is(typeof(a) == uint));
}
