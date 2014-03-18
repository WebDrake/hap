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

import std.range, std.stdio, std.traits;

/**
 * Reads randomness from an infinite filestream.  In practice this
 * will typically be used to read from system sources of randomness
 * such as (on Posix) $(D /dev/random) or $(D /dev/urandom).
 *
 * The random numbers generated will be unsigned integers of type
 * $(D T) in the range [$(D T.min), $(D T.max)].  A good source of
 * randomness should ensure these are uniformly distributed.
 *
 * It is the responsibility of the user to ensure that the specified
 * source of randomness indeed contains sufficient data to serve the
 * requirements of their program.
 */
final class RandomFileStream(string filename, T)
    if (isIntegral!T && isUnsigned!T)
{
  private:
    File _dev;
    T _value;

  public:
    alias source = filename;
    enum T min = T.min;
    enum T max = T.max;

    enum bool isUniformRandom = true;

    static assert (this.min == 0);

    this()
    {
        _dev = File(filename, "r");
        popFront();
    }

    /// Range primitives
    enum bool empty = false;

    /// ditto
    T front() @property @safe const nothrow pure
    {
        return _value;
    }

    /// ditto
    void popFront()
    {
        _dev.rawRead((&_value)[0 .. 1]);
    }
}

version (Posix)
{
    /**
     * Generates uniformly distributed random numbers in the interval
     * [$(D T.min), $(D T.max)] using $(D /dev/random) as the source
     * of randomness.
     *
     * Caution should be taken when using this as $(D /dev/random) is
     * a blocking device.
     */
    template DevRandom(T)
        if (isIntegral!T && isUnsigned!T)
    {
        alias DevRandom = RandomFileStream!("/dev/random", T);
    }

    /**
     * Generates uniformly distributed random numbers in the interval
     * [$(D T.min), $(D T.max)] using $(D /dev/urandom) as the source
     * of randomness.
     */
    template DevURandom(T)
        if (isIntegral!T && isUnsigned!T)
    {
        alias DevURandom = RandomFileStream!("/dev/urandom", T);
    }
}

unittest
{
    import std.stdio, std.typetuple;
    writeln("Imported std.random2.device!");

    foreach (Dev; TypeTuple!(DevRandom, DevURandom))
    {
        foreach (T; TypeTuple!(ushort, uint, ulong))
        {
            static assert(isUniformRNG!(Dev!T));

            auto rnd = new Dev!T;

            writeln("Reading random numbers of type " ~ T.stringof ~ " from " ~ rnd.source);

            foreach (var; rnd.take(20))
            {
                writeln(var);
            }
            writeln;

            writeln("Generating U[0, 1) from " ~ rnd.source ~ " variates of type " ~ T.stringof);

            import std.random2.distribution;

            auto udist = uniformDistribution(0.0, 1.0, rnd);

            foreach (u; udist.take(20))
            {
                writeln(u);
            }
            writeln;
        }
    }
}
