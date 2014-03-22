// Written in the D programming language.

/**
 * Implements algorithms for generating random numbers drawn from
 * different statistical distributions.  Where possible, each random
 * distribution is provided in two different forms:
 *
 * $(UL
 *   $(LI as a function, which takes as input the distribution
 *        parameters and a uniform RNG, and returns a single value
 *        drawn from the distribution, and)
 *   $(LI as a range object, which wraps a uniform RNG instance and
 *        transforms its output into numbers drawn from the specified
 *        distribution.)
 * )
 *
 * Typical reasons for rejecting a function implementation include
 * the function needing to hold state between calls to achieve
 * adequate performance, or the function needing to allocate memory
 * with each call.
 *
 * As with random number generators, the random distribution range
 * objects implemented here are final classes in order to ensure
 * reference semantics.  They also assume reference type semantics on
 * the part of the RNGs that they wrap: user-supplied value-type RNGs
 * may produce unexpected and incorrect behaviour when combined with
 * these objects.
 *
 * Copyright: © 2014 Joseph Rushton Wakeling
 *
 * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0).
 *
 * Authors: $(WEB erdani.org, Andrei Alexandrescu),
 *          Chris Cain,
 *          $(WEB braingam.es, Joseph Rushton Wakeling)
 *
 * Source: $(PHOBOSSRC std/random2/_distribution.d)
 */
module std.random2.distribution;

import std.random2.generator, std.random2.traits;

import std.range, std.traits;

// dice
/**
 * Rolls a random die with relative probabilities stored in $(D proportions).
 * Returns the index in $(D proportions) that was chosen.
 *
 * Example:
 * ----
 * auto x = dice(0.5, 0.5);   // x is 0 or 1 in equal proportions
 * auto y = dice(50, 50);     // y is 0 or 1 in equal proportions
 * auto z = dice(70, 20, 10); // z is 0 70% of the time, 1 20% of the time,
 *                            // and 2 10% of the time
 * ----
 *
 * The range counterpart of $(D dice) is the $(D DiscreteDistribution) class.
 */
size_t dice(UniformRNG, Num)(ref UniformRNG rng, Num[] proportions...)
    if (isNumeric!Num && isForwardRange!UniformRNG)
{
    return diceImpl(rng, proportions);
}

/// ditto
size_t dice(UniformRNG, Range)(ref UniformRNG rng, Range proportions)
    if (isUniformRNG!UniformRNG && isForwardRange!Range &&
        isNumeric!(ElementType!Range) && !isArray!Range)
{
    return diceImpl(rng, proportions);
}

/// ditto
size_t dice(Range)(Range proportions)
    if (isForwardRange!Range && isNumeric!(ElementType!Range) && !isArray!Range)
{
    return diceImpl(rndGen, proportions);
}

/// ditto
size_t dice(Num)(Num[] proportions...)
    if (isNumeric!Num)
{
    return diceImpl(rndGen, proportions);
}

private size_t diceImpl(UniformRNG, Range)(ref UniformRNG rng, Range proportions)
    if (isUniformRNG!UniformRNG && isForwardRange!Range && isNumeric!(ElementType!Range))
{
    import std.algorithm, std.exception, std.random2.distribution;
    double sum = reduce!("(assert(b >= 0), a + b)")(0.0, proportions.save);
    enforce(sum > 0, "Proportions in a dice cannot sum to zero");
    immutable point = uniform(0.0, sum, rng);
    assert(point < sum);
    auto mass = 0.0;

    size_t i = 0;
    foreach (e; proportions)
    {
        mass += e;
        if (point < mass)
        {
            return i;
        }
        i++;
    }
    // this point should not be reached
    assert(false);
}

unittest
{
    auto gen = new Random(unpredictableSeed);
    auto i = dice(gen, 0.0, 100.0);
    assert(i == 1);
    i = dice(gen, 100.0, 0.0);
    assert(i == 0);

    i = dice(100U, 0U);
    assert(i == 0);
}

/**
 * The range equivalent of $(D dice), each element of which is the
 * result of the roll of a random die with relative probabilities
 * stored in the range $(D proportions).  Each successive value of
 * $(D front) reflects the index in $(D proportions) that was chosen.
 *
 * This offers a superior option in the event of making many rolls
 * of the same die, as the sum and CDF of $(D proportions) only needs
 * to be calculated upon construction and not with each call.
 */
final class DiscreteDistribution(SearchPolicy search, T, UniformRNG)
{
  private:
    SortedRange!(immutable(T)[]) _cumulativeProportions;
    UniformRNG _rng;
    size_t _value;

  public:
    enum bool isRandomDistribution = true;

    this(Range)(UniformRNG rng, Range proportions)
        if (isInputRange!Range && isNumeric!(ElementType!Range))
    in
    {
        assert(!proportions.empty, "Proportions of discrete distribution cannot be empty.");
    }
    body
    {
        import std.exception;

        _rng = rng;

        static if (isImplicitlyConvertible!(typeof(proportions.array), T[]))
        {
            T[] prop = proportions.array;
        }
        else
        {
            import std.conv;
            T[] prop = proportions.map!(to!T).array;
        }

        alias A = Select!(isFloatingPoint!T, real, T);

        A accumulator = 0;

        foreach(ref e; prop)
        {
            assert(e >= 0, "Proportions of discrete distribution cannot be negative.");

            debug
            {
                A preAccumulation = accumulator;
            }

            accumulator += e;

            debug
            {
                static if (isIntegral!T)
                {
                    enforce(preAccumulation <= accumulator, "Integer overflow detected!");
                }
                else static if (isFloatingPoint!T)
                {
                    if (e > 0)
                    {
                        enforce(accumulator > preAccumulation, "Floating point rounding error detected!");
                    }
                }
                else
                {
                    static assert(0);
                }
            }

            e = accumulator;
        }

        _cumulativeProportions = assumeSorted(assumeUnique(prop)[]);

        popFront();
    }

    this(typeof(this) that)
    {
        this._cumulativeProportions = that._cumulativeProportions.save;
        this._rng = that._rng;
        this._value = that._value;
    }

    /// Range primitives
    enum bool empty = false;

    /// ditto
    size_t front() @property @safe const nothrow pure
    {
        return _value;
    }

    /// ditto
    void popFront()
    {
        immutable sum = _cumulativeProportions.back;
        immutable point = uniform!"[)"(0, sum, _rng);
        assert(point < sum);
        _value = _cumulativeProportions.length - _cumulativeProportions.upperBound!search(point).length;
    }

    /// ditto
    static if (isForwardRange!UniformRNG)
    {
        typeof(this) save() @property
        {
            auto ret = new typeof(this)(this);
            ret._rng = this._rng.save;
            return ret;
        }
    }
}

/// ditto
auto discreteDistribution(SearchPolicy search = SearchPolicy.binarySearch, UniformRNG, Range)
                         (UniformRNG rng, Range proportions)
    if (isInputRange!Range && isNumeric!(ElementType!Range) && isUniformRNG!UniformRNG)
{
    alias E = Unqual!(ElementType!Range);

    static if (isIntegral!E)
    {
        static if (is(Largest!(ushort, Unsigned!E) == ushort))
        {
            alias T = uint;
        }
        else static if (is(Largest!(ulong, Unsigned!E) == ulong))
        {
            alias T = ulong;
        }
        else
        {
            static assert(false);
        }
    }
    else static if (isFloatingPoint!E)
    {
        alias T = Largest!(double, E);
    }

    return new DiscreteDistribution!(search, T, UniformRNG)(rng, proportions);
}

/// ditto
auto discreteDistribution(SearchPolicy search = SearchPolicy.binarySearch, Range)
                         (Range proportions)
    if (isInputRange!Range && isNumeric!(ElementType!Range))
{
    return discreteDistribution(rndGen, proportions);
}

/// ditto
auto discreteDistribution(SearchPolicy search = SearchPolicy.binarySearch, UniformRNG, Num)
                         (UniformRNG rng, Num[] proportions...)
    if (isUniformRNG!UniformRNG && isNumeric!Num)
{
    return discreteDistribution(rng, proportions);
}

/// ditto
auto discreteDistribution(SearchPolicy search = SearchPolicy.binarySearch, Num)
                         (Num[] proportions...)
    if (isNumeric!Num)
{
    return discreteDistribution(rndGen, proportions);
}

unittest
{
    import std.random2.device, std.stdio;

    foreach (UniformRNG; TypeTuple!(UniformRNGTypes, DevURandom!ulong))
    {
        auto ddist = discreteDistribution(25, 50, 25);
        size_t[3] prop;

        prop[] = 0;

        writeln("Discrete Distribution (", UniformRNG.stringof, "): 25:50:25");
        foreach (d; ddist.take(100_000))
        {
            prop[d]++;
        }
        writeln(prop);
    }
}

/**
 * Provides an infinite range of random numbers distributed according to the
 * normal (Gaussian) distribution with mean $(D mu) and standard deviation
 * $(D sigma).  The version that does not receive a specified random number
 * generator uses the default generator $(D rndGen) as its source of
 * randomness.
 */
final class NormalDistribution(T, UniformRNG)
    if (isFloatingPoint!T && isUniformRNG!UniformRNG)
{
  private:
    alias NormalEngine = NormalEngineBoxMuller;
    NormalEngine!T _engine;
    UniformRNG _rng;
    T _value;

  public:
    enum bool isRandomDistribution = true;
    immutable T mean;
    immutable T stdev;

    this(T mu, T sigma, UniformRNG rng)
    {
        import std.exception;
        enforce(sigma > 0);
        mean = mu;
        stdev = sigma;
        _rng = rng;
        popFront();
    }

    this(typeof(this) that)
    {
        this(that.mean, that.stdev, that._rng);
    }

    /// Range primitives.
    enum bool empty = false;

    /// ditto
    T front() @property @safe const nothrow pure
    {
        return _value;
    }

    /// ditto
    void popFront()
    {
        _value = _engine(this.mean, this.stdev, _rng);
    }

    /// ditto
    static if (isForwardRange!UniformRNG)
    {
        typeof(this) save() @property
        {
            auto ret = new typeof(this)(this);
            ret._rng = this._rng.save;
            return ret;
        }
    }
}

/// ditto
auto normalDistribution(T1, T2, UniformRNG)
                       (T1 mu, T2 sigma, UniformRNG rng)
    if (isNumeric!T1 && isNumeric!T2 && isUniformRNG!UniformRNG)
{
    static if (isFloatingPoint!(CommonType!(T1, T2)))
    {
        alias T = CommonType!(T1, T2);
    }
    else
    {
        alias T = double;
    }

    return new NormalDistribution!(T, UniformRNG)(mu, sigma, rng);
}

/// ditto
auto normalDistribution(T1, T2)
                       (T1 mu, T2 sigma)
    if (isNumeric!T1 && isNumeric!T2)
{
    return normalDistribution!(T1, T2, Random)(mu, sigma, rndGen);
}

unittest
{
    // check type rules for NormalDistribution
    {
        auto ndist = normalDistribution(0, 1);
        static assert(is(typeof(ndist.front) == double));
    }
    {
        auto ndist = normalDistribution(0.0f, 1);
        static assert(is(typeof(ndist.front) == float));
    }
    {
        auto ndist = normalDistribution(0.0, 1);
        static assert(is(typeof(ndist.front) == double));
    }
    {
        auto ndist = normalDistribution(0.0L, 1);
        static assert(is(typeof(ndist.front) == real));
    }

    // check save works effectively
    {
        auto ndist = normalDistribution(3.29, 7.64);
        auto ndist2 = ndist.save;
        assert(ndist2 !is ndist);

        ndist.popFrontN(10);
        assert(ndist2.front != ndist.front);
        ndist2.popFrontN(10);
        assert(ndist2.front == ndist.front);
    }
}

/**
 * Generates random numbers drawn from a normal (Gaussian) distribution, using
 * the Box-Muller Transform method.
 *
 * This implementation of Box-Muller closely follows that of its counterpart
 * in the Boost.Random C++ library and should produce matching results aside
 * from discrepancies that arise out of differences in floating-point precision.
 */
private struct NormalEngineBoxMuller(T)
    if (isFloatingPoint!T)
{
  private:
    bool _valid = false;
    T _rho, _r1, _r2;

  public:
    /**
     * Generates a single random number drawn from a normal distribution with
     * mean $(D mu) and standard deviation $(D sigma), using $(D rng) as the
     * source of randomness.
     */
    T opCall(UniformRNG)(in T mu, in T sigma, ref UniformRNG rng)
        if (isUniformRNG!UniformRNG)
    {
        import std.math;

        _valid = !_valid;

        if (_valid)
        {
            /* N.B. Traditional Box-Muller asks for random numbers
             * in (0, 1], which uniform() can readily supply.  We
             * instead generate numbers in [0, 1) and use 1 - num
             * to match the output of Boost.Random.
             */
            _r1 = uniform!("[)", T, T, UniformRNG)(0, 1, rng);
            _r2 = uniform!("[)", T, T, UniformRNG)(0, 1, rng);
            _rho = sqrt(-2 * log(1 - _r2));

            return _rho * cos(2 * PI * _r1) * sigma + mu;
        }
        else
        {
            return _rho * sin(2 * PI * _r1) * sigma + mu;
        }
    }
}

unittest
{
    NormalEngineBoxMuller!double engine;
    auto rng1 = new Random(unpredictableSeed);
    auto rng2 = rng1.save;
    auto rng3 = rng1.save;
    double mu = 6.5, sigma = 3.2;

    /* The Box-Muller engine produces variates a pair at
     * a time.  We verify this is true by using a pair of
     * pseudo-random number generators sharing the same
     * initial state.
     */
    auto a1 = engine(mu, sigma, rng1);
    auto b2 = engine(mu, sigma, rng2);

    // verify that 1st RNG has been called but 2nd has not
    assert(rng3.front != rng1.front);
    assert(rng3.front == rng2.front);

    /* Now, calling with the RNG order reversed should
     * produce the same results: only rng2 will get called
     * this time.
     */
    auto a2 = engine(mu, sigma, rng2);
    auto b1 = engine(mu, sigma, rng1);

    assert(a1 == a2);
    assert(b1 == b2);
    assert(rng2.front == rng1.front);
    assert(rng3.front != rng2.front);

    // verify that the RNGs have each been called twice
    rng3.popFrontN(2);
    assert(rng3.front == rng2.front);
}

/**
 * Generates a number between $(D a) and $(D b). The $(D boundaries)
 * parameter controls the shape of the interval (open vs. closed on
 * either side). Valid values for $(D boundaries) are $(D "[]"), $(D
 * "$(LPAREN)]"), $(D "[$(RPAREN)"), and $(D "()"). The default interval
 * is closed to the left and open to the right. The version that does not
 * take $(D rng) uses the default generator $(D rndGen).
 *
 * Example:
 *
 * ----
 * auto gen = Random(unpredictableSeed);
 * // Generate an integer in [0, 1023]
 * auto a = uniform(0, 1024, gen);
 * // Generate a float in [0, 1$(RPAREN)
 * auto a = uniform(0.0f, 1.0f, gen);
 * ----
 */
auto uniform(string boundaries = "[)", T1, T2)
            (T1 a, T2 b)
    if (!is(CommonType!(T1, T2) == void))
{
    return uniform!(boundaries, T1, T2, Random)(a, b, rndGen);
}

unittest
{
    //MinstdRand0 gen;
    auto gen = new Mt19937;
    foreach (i; 0 .. 20)
    {
        auto x = uniform(0.0, 15.0, gen);
        assert(0 <= x && x < 15);
    }
    foreach (i; 0 .. 20)
    {
        auto x = uniform!"[]"('a', 'z', gen);
        assert('a' <= x && x <= 'z');
    }

    foreach (i; 0 .. 20)
    {
        auto x = uniform('a', 'z', gen);
        assert('a' <= x && x < 'z');
    }

    foreach(i; 0 .. 20)
    {
        immutable ubyte a = 0;
            immutable ubyte b = 15;
        auto x = uniform(a, b, gen);
            assert(a <= x && x < b);
    }
}

// Implementation of uniform for floating-point types
/// ditto
auto uniform(string boundaries = "[)", T1, T2, UniformRNG)
            (T1 a, T2 b, ref UniformRNG rng)
    if (isFloatingPoint!(CommonType!(T1, T2)) && isUniformRNG!UniformRNG)
out (result)
{
    static if (boundaries[0] == '(')
    {
        assert(a < result);
    }
    else static if (boundaries[0] == '[')
    {
        assert(a <= result);
    }

    static if (boundaries[1] == ')')
    {
        assert(result < b);
    }
    else static if (boundaries[1] == ']')
    {
        assert(result <= b);
    }
}
body
{
    import std.exception, std.math, std.string : format;
    alias Unqual!(CommonType!(T1, T2)) NumberType;
    static if (boundaries[0] == '(')
    {
        NumberType _a = nextafter(cast(NumberType) a, NumberType.infinity);
    }
    else
    {
        NumberType _a = a;
    }
    static if (boundaries[1] == ')')
    {
        NumberType _b = nextafter(cast(NumberType) b, -NumberType.infinity);
    }
    else
    {
        NumberType _b = b;
    }
    enforce(_a <= _b,
            format("std.random.uniform(): invalid bounding interval %s%s, %s%s",
                   boundaries[0], a, b, boundaries[1]));
    NumberType result =
        _a + (_b - _a) * cast(NumberType) (rng.front - rng.min)
        / (rng.max - rng.min);
    rng.popFront();
    return result;
}

// Implementation of uniform for integral types
/// ditto
auto uniform(string boundaries = "[)", T1, T2, UniformRNG)
            (T1 a, T2 b, ref UniformRNG rng)
    if ((isIntegral!(CommonType!(T1, T2)) || isSomeChar!(CommonType!(T1, T2)))
        && isUniformRNG!UniformRNG)
out (result)
{
    static if (boundaries[0] == '(')
    {
        assert(a < result);
    }
    else static if (boundaries[0] == '[')
    {
        assert(a <= result);
    }

    static if (boundaries[1] == ')')
    {
        assert(result < b);
    }
    else static if (boundaries[1] == ']')
    {
        assert(result <= b);
    }
}
body
{
    import std.conv, std.exception, std.string : format;
    alias Unqual!(CommonType!(T1, T2)) ResultType;
    // We handle the case "[)' as the common case, and we adjust all
    // other cases to fit it.
    static if (boundaries[0] == '(')
    {
        enforce(cast(ResultType) a < ResultType.max,
                format("std.random.uniform(): invalid left bound %s", a));
        ResultType min = cast(ResultType) a + 1;
    }
    else
    {
        ResultType min = a;
    }
    static if (boundaries[1] == ']')
    {
        enforce(min <= cast(ResultType) b,
                format("std.random.uniform(): invalid bounding interval %s%s, %s%s",
                       boundaries[0], a, b, boundaries[1]));
        if (b == ResultType.max && min == ResultType.min)
        {
            // Special case - all bits are occupied
            return .uniform!ResultType(rng);
        }
        auto count = unsigned(b - min) + 1u;
        static assert(count.min == 0);
    }
    else
    {
        enforce(min < cast(ResultType) b,
                format("std.random.uniform(): invalid bounding interval %s%s, %s%s",
                       boundaries[0], a, b, boundaries[1]));
        auto count = unsigned(b - min);
        static assert(count.min == 0);
    }
    assert(count != 0);
    if (count == 1) return min;
    alias typeof(count) CountType;
    static assert(CountType.min == 0);
    auto bucketSize = 1u + (CountType.max - count + 1) / count;
    CountType r;
    do
    {
        r = cast(CountType) (uniform!CountType(rng) / bucketSize);
    }
    while (r >= count);
    return cast(typeof(return)) (min + r);
}

unittest
{
    import std.conv, std.typetuple;
    auto gen = new Mt19937(unpredictableSeed);
    static assert(isForwardRange!(typeof(gen)));

    auto a = uniform(0, 1024, gen);
    assert(0 <= a && a <= 1024);
    auto b = uniform(0.0f, 1.0f, gen);
    assert(0 <= b && b < 1, to!string(b));
    auto c = uniform(0.0, 1.0);
    assert(0 <= c && c < 1);

    foreach(T; TypeTuple!(char, wchar, dchar, byte, ubyte, short, ushort,
                          int, uint, long, ulong, float, double, real))
    {
        T lo = 0, hi = 100;
        T init = uniform(lo, hi);
        size_t i = 50;
        while (--i && uniform(lo, hi) == init) {}
        assert(i > 0);
    }
}

/**
 * Generates a uniformly-distributed number in the range $(D [T.min, T.max])
 * for any integral type $(D T). If no random number generator is passed,
 * uses the default $(D rndGen).
 */
auto uniform(T, UniformRNG)(ref UniformRNG rng)
    if (!is(T == enum) && (isIntegral!T || isSomeChar!T)
        && isUniformRNG!UniformRNG)
{
    auto r = rng.front;
    rng.popFront();
    static if (T.sizeof <= r.sizeof)
    {
        return cast(T) r;
    }
    else
    {
        static assert(T.sizeof == 8 && r.sizeof == 4);
        T r1 = rng.front | (cast(T)r << 32);
        rng.popFront();
        return r1;
    }
}

/// ditto
auto uniform(T)()
    if (!is(T == enum) && (isIntegral!T || isSomeChar!T))
{
    return uniform!T(rndGen);
}

unittest
{
    import std.typetuple;
    foreach(T; TypeTuple!(char, wchar, dchar, byte, ubyte, short, ushort,
                          int, uint, long, ulong))
    {
        T init = uniform!T();
        size_t i = 50;
        while (--i && uniform!T() == init) {}
        assert(i > 0);
        assert(i < 50);
    }
}

/**
 * Provides an infinite sequence of random numbers uniformly distributed between
 * $(D a) and $(D b).  The $(D boundaries) parameter controls the shape of the
 * interval (open vs. closed on either side).  Valid values for $(D boundaries)
 * are $(D "[]"), $(D "$(LPAREN)]"), $(D "[$(RPAREN)"), and $(D "()"). The
 * default interval is closed to the left and open to the right. The version
 * that does not receive a specified random number generator uses the default
 * generator $(D rndGen).
 */
final class UniformDistribution(string boundaries = "[)", T, UniformRNG)
    if (isNumeric!T && isUniformRNG!UniformRNG)
{
  private:
    UniformRNG _rng;
    T _value;

  public:
    enum bool isRandomDistribution = true;

    immutable T min;
    immutable T max;

    this(T a, T b, UniformRNG rng)
    {
        import std.exception;
        enforce(a < b);
        min = a;
        max = b;
        _rng = rng;
        popFront();
    }

    this(typeof(this) that)
    {
        this(that.min, that.max, that._rng);
    }

    /// Range primitives.
    enum bool empty = false;

    /// ditto
    T front() @property @safe const nothrow pure
    {
        return _value;
    }

    /// ditto
    void popFront()
    {
        _value = uniform!(boundaries, T, T, UniformRNG)(min, max, _rng);
    }

    /// ditto
    static if (isForwardRange!UniformRNG)
    {
        typeof(this) save() @property
        {
            auto ret = new typeof(this)(this);
            ret._rng = this._rng.save;
            return ret;
        }
    }
}

/// ditto
auto uniformDistribution(string boundaries = "[)", T1, T2, UniformRNG)
                        (T1 a, T2 b, ref UniformRNG rng)
    if (isNumeric!(CommonType!(T1, T2)) && isUniformRNG!UniformRNG)
{
    alias T = CommonType!(T1, T2);
    return new UniformDistribution!(boundaries, T, UniformRNG)(a, b, rng);
}

/// ditto
auto uniformDistribution(string boundaries = "[)", T1, T2)
                        (T1 a, T2 b)
    if (isNumeric!(CommonType!(T1, T2)))
{
    alias T = CommonType!(T1, T2);
    return new UniformDistribution!(boundaries, T, Random)(a, b, rndGen);
}

unittest
{
    static assert(isRandomDistribution!(UniformDistribution!("[)", int, Random)));
    static assert(isRandomDistribution!(UniformDistribution!("(]", double, Random)));

    // check distribution boundaries function OK for floating-point
    {
        auto udist = uniformDistribution(3.2, 5.9);
        foreach (u; udist.take(100))
        {
            assert(udist.min <= u);
            assert(u < udist.max);
        }
    }
    {
        auto udist = uniformDistribution!"()"(5.7, 8.9);
        foreach (u; udist.take(100))
        {
            assert(udist.min < u);
            assert(u < udist.max);
        }
    }
    {
        auto udist = uniformDistribution!"(]"(3.7L, 9.2L);
        foreach (u; udist.take(100))
        {
            assert(udist.min < u);
            assert(u <= udist.max);
        }
    }
    {
        auto udist = uniformDistribution!"[]"(2.1, 9.9);
        foreach (u; udist.take(100))
        {
            assert(udist.min <= u);
            assert(u <= udist.max);
        }
    }

    // check distribution boundaries function OK for integers
    {
        auto udist = uniformDistribution(0, 10);
        size_t eqMin = 0;
        foreach (u; udist.take(100))
        {
            assert(udist.min <= u);
            assert(u < udist.max);
            if (u == udist.min) ++eqMin;
        }
        assert(eqMin > 0);
    }
    {
        auto udist = uniformDistribution!"()"(0, 10);
        foreach (u; udist.take(100))
        {
            assert(udist.min < u);
            assert(u < udist.max);
        }
    }
    {
        auto udist = uniformDistribution!"(]"(0, 10);
        size_t eqMax = 0;
        foreach (u; udist.take(100))
        {
            assert(udist.min < u);
            assert(u <= udist.max);
            if (u == udist.max) ++eqMax;
        }
        assert(eqMax > 0);
    }
    {
        auto udist = uniformDistribution!"[]"(0, 10);
        size_t eqMin = 0, eqMax = 0;
        foreach (u; udist.take(100))
        {
            assert(udist.min <= u);
            assert(u <= udist.max);
            if (u == udist.min) ++eqMin;
            if (u == udist.max) ++eqMax;
        }
        assert(eqMin > 0);
        assert(eqMax > 0);
    }

    // check that save works properly
    {
        auto udist = uniformDistribution(1.1, 3.3);
        auto udist2 = udist.save;
        assert(udist2 !is udist);
        udist.popFrontN(20);
        assert(udist2.front != udist.front);
        udist2.popFrontN(20);
        assert(udist2.front == udist.front);
    }
}

/**
 * Generates a uniformly-distributed floating point number
 * of type $(D T) in the range [0, 1).  If no random number
 * generator is specified, the default RNG $(D rndGen) will
 * be used as the source of randomness.
 *
 * uniform01 offers a faster generation of random variates
 * than the equivalent uniform(0.0, 1.0) and so may be
 * preferred for some applications.
 */
T uniform01(T = double)()
    if (isFloatingPoint!T)
{
    return uniform01!T(rndGen);
}

/// ditto
T uniform01(T = double, UniformRNG)(UniformRNG rng)
    if (isFloatingPoint!T)
out (result)
{
    assert(result < 1);
}
body
{
    import std.string : format;
    alias R = typeof(rng.front);
    static if (isIntegral!R)
    {
        enum T denom = 1 / ((cast(T) 1) + rng.max - rng.min);
    }
    else static if (isFloatingPoint!R)
    {
        enum T denom = 1 / ((cast(T) 1) + rng.max - rng.min);
    }
    else
    {
        static assert(false);
    }

    T u = (rng.front - rng.min) * denom;

    rng.popFront();

    return u;
}

unittest
{
    import std.typetuple;
    foreach (UniformRNG; UniformRNGTypes)
    {

        foreach (T; TypeTuple!(float, double, real))
        {
            UniformRNG rng = new UniformRNG(unpredictableSeed);

            auto a = uniform01();
            assert(is(typeof(a) == double));
            assert(0 <= a && a < 1);

            auto b = uniform01(rng);
            assert(is(typeof(a) == double));
            assert(0 <= b && b < 1);

            auto c = uniform01!T();
            assert(is(typeof(c) == T));
            assert(0 <= c && c < 1);

            auto d = uniform01!T(rng);
            assert(is(typeof(d) == T));
            assert(0 <= d && d < 1);

            T init = uniform01!T(rng);
            size_t i = 50;
            while (--i && uniform01!T(rng) == init) {}
            assert(i > 0);
            assert(i < 50);
        }
    }
}

/**
 * Provides an infinite sequence of random numbers uniformly distributed in the
 * interval [0, 1).  If no RNG is specified, $(D uniformDistribution) will use
 * the default generator $(D rndGen).
 */
final class Uniform01Distribution(T, UniformRNG)
    if (isFloatingPoint!T && isUniformRNG!UniformRNG)
{
  private:
    UniformRNG _rng;
    T _value;

  public:
    enum T min = 0;
    enum T max = 1;
    enum bool isRandomDistribution = true;
    enum bool isUniformRandom = true;

    this(UniformRNG rng)
    {
        _rng = rng;
        popFront();
    }

    this(typeof(this) that)
    {
        this(that._rng);
    }

    /// Range primitives.
    enum bool empty = false;

    /// ditto
    T front() @property @safe const nothrow pure
    {
        return _value;
    }

    /// ditto
    void popFront()
    {
        _value = uniform01!(T, UniformRNG)(_rng);
    }

    /// ditto
    static if (isForwardRange!UniformRNG)
    {
        typeof(this) save() @property
        {
            auto ret = new typeof(this)(this);
            ret._rng = this._rng.save;
            return ret;
        }
    }
}

/// ditto
auto uniform01Distribution(T = double, UniformRNG)(ref UniformRNG rng)
    if (isFloatingPoint!T && isUniformRNG!UniformRNG)
{
    return new Uniform01Distribution!(T, UniformRNG)(rng);
}

/// ditto
auto uniform01Distribution(T = double)()
    if (isFloatingPoint!T)
{
    return new Uniform01Distribution!(T, Random)(rndGen);
}

unittest
{
    import std.stdio;
    auto u01 = uniform01Distribution();

    foreach (immutable _; 0 .. 1_000_000)
    {
        auto u = uniform01(u01);
    }

    auto u0101 = uniform01Distribution(u01);

    foreach (u; u0101.take(1_000_000))
    {
    }

    import std.random2.adaptor;

    iota(100).sample(5, u0101).writeln;
}
