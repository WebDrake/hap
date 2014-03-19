import std.datetime, std.random, std.range, std.stdio, std.typetuple;

// $(D TypeTuple) of all uniform RNGs defined in std.random.
alias UniformRNGTypes =
    TypeTuple!(MinstdRand0, MinstdRand,
               Mt19937,
               Xorshift32, Xorshift64, Xorshift128, Xorshift160, Xorshift192);

void main()
{
    enum k = 10_000;
    enum l = 1_000_000;
    enum m = 10_000_000;
    enum n = 100_000_000;

    foreach (UniformRNG; UniformRNGTypes)
    {
        StopWatch watch;
        auto rng = UniformRNG(unpredictableSeed);
        writeln("Testing random generator: ", typeof(rng).stringof);
        writeln("\tfront: ", rng.front);
        watch.start();
        rng.popFrontN(n);
        watch.stop();
        writeln("\tfront after ", n, " pops: ", rng.front);
        writeln("\tcalculated in ", watch.peek.msecs, " ms.");

        auto proportions = cycle([1, 2, 3, 4, 5]).take(100).array;

        watch.reset();
        watch.start();
        foreach (immutable _; 0 .. m)
        {
            auto d = dice(rng, proportions);
        }
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to perform ", m, " calls to dice.");

        watch.reset();
        watch.start();
        foreach (immutable _; 0 .. m)
        {
            auto u = uniform(0.0, 1.0, rng);
        }
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to perform ", m, " calls to uniform.");

        watch.reset();
        watch.start();
        foreach (immutable _; 0 .. k)
        {
            auto cover = randomCover(proportions, rng);
            foreach (c; cover)
            {
                auto val = c;
            }
        }
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to generate ", k, " random covers of ",
                proportions.length, " elements.");

        watch.reset();
        watch.start();
        foreach (immutable _; 0 .. l)
        {
            auto sample = randomSample(proportions, 5, rng);
            foreach (s; sample)
            {
                auto val = s;
            }
        }
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to generate ", l, " random samples of 5 out of ",
                proportions.length, " elements.");

        watch.reset();
        watch.start();
        foreach (immutable _; 0 .. l)
        {
            randomShuffle(proportions, rng);
        }
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to generate ", l, " in-place random shuffles of ",
                proportions.length, " elements.");

        writeln();
    }
}
