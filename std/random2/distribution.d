module std.random2.distribution;

import std.random2.generator, std.random2.traits;

import std.range, std.traits;

/**
Generates a number between $(D a) and $(D b). The $(D boundaries)
parameter controls the shape of the interval (open vs. closed on
either side). Valid values for $(D boundaries) are $(D "[]"), $(D
"$(LPAREN)]"), $(D "[$(RPAREN)"), and $(D "()"). The default interval
is closed to the left and open to the right. The version that does not
take $(D urng) uses the default generator $(D rndGen).

Example:

----
auto gen = Random(unpredictableSeed);
// Generate an integer in [0, 1023]
auto a = uniform(0, 1024, gen);
// Generate a float in [0, 1$(RPAREN)
auto a = uniform(0.0f, 1.0f, gen);
----
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
            (T1 a, T2 b, ref UniformRNG urng)
    if (isFloatingPoint!(CommonType!(T1, T2)) && isUniformRNG!UniformRNG)
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
        _a + (_b - _a) * cast(NumberType) (urng.front - urng.min)
        / (urng.max - urng.min);
    urng.popFront();
    return result;
}

// Implementation of uniform for integral types
auto uniform(string boundaries = "[)", T1, T2, UniformRNG)
            (T1 a, T2 b, ref UniformRNG urng)
    if ((isIntegral!(CommonType!(T1, T2)) || isSomeChar!(CommonType!(T1, T2)))
        && isUniformRNG!UniformRNG)
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
            return .uniform!ResultType(urng);
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
        r = cast(CountType) (uniform!CountType(urng) / bucketSize);
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
Generates a uniformly-distributed number in the range $(D [T.min,
T.max]) for any integral type $(D T). If no random number generator is
passed, uses the default $(D rndGen).
 */
auto uniform(T, UniformRNG)(ref UniformRNG urng)
    if (!is(T == enum) && (isIntegral!T || isSomeChar!T)
        && isUniformRNG!UniformRNG)
{
    auto r = urng.front;
    urng.popFront();
    static if (T.sizeof <= r.sizeof)
    {
        return cast(T) r;
    }
    else
    {
        static assert(T.sizeof == 8 && r.sizeof == 4);
        T r1 = urng.front | (cast(T)r << 32);
        urng.popFront();
        return r1;
    }
}

/// Ditto
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

auto uniformDistribution(string boundaries = "[)", T1, T2, UniformRNG)
                        (T1 a, T2 b, ref UniformRNG rng)
    if (isNumeric!(CommonType!(T1, T2)) && isUniformRNG!UniformRNG)
{
    alias T = CommonType!(T1, T2);
    return new UniformDistribution!(boundaries, T, UniformRNG)(a, b, rng);
}

auto uniformDistribution(string boundaries = "[)", T1, T2)
                        (T1 a, T2 b)
    if (isNumeric!(CommonType!(T1, T2)))
{
    alias T = CommonType!(T1, T2);
    return new UniformDistribution!(boundaries, T, Random)(a, b, rndGen);
}

unittest
{
    import std.stdio;

    static assert(isRandomDistribution!(UniformDistribution!("[)", int, Random)));
    static assert(isRandomDistribution!(UniformDistribution!("(]", double, Random)));

    auto udist = uniformDistribution(3.2, 5.9);

    foreach (u; udist.take(100))
    {
        writeln("U(", udist.min, ", ", udist.max, ") = ", u);
        assert(udist.min <= u);
        assert(u <= udist.max);
    }

    auto udist2 = udist.save;
    assert(udist2 !is udist);
    udist.popFrontN(20);
    assert(udist2.front != udist.front);
    udist2.popFrontN(20);
    assert(udist2.front == udist.front);
}
