import std.datetime, std.random2, std.range, std.stdio;

void main()
{
    enum k = 10_000;
    enum l = 1_000_000;
    enum m = 10_000_000;
    enum n = 100_000_000;

    foreach (UniformRNG; UniformRNGTypes)
    {
        StopWatch watch;
        auto rng = new UniformRNG(unpredictableSeed);
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

        auto ddist = discreteDistribution(rng, proportions);

        watch.reset();
        watch.start();
        ddist.popFrontN(m);
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to perform ", m, " pops of ", typeof(ddist).stringof);

        auto ndist = normalDistribution(0.0, 1.0, rng);

        watch.reset();
        watch.start();
        ndist.popFrontN(m);
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to perform ", m, " pops of ",
                typeof(ndist).stringof);

        watch.reset();
        watch.start();
        foreach (immutable _; 0 .. m)
        {
            auto u = uniform(0.0, 1.0, rng);
        }
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to perform ", m, " calls to uniform.");

        auto udist = uniformDistribution(0.0, 1.0, rng);

        watch.reset();
        watch.start();
        udist.popFrontN(m);
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to perform ", m, " pops of ", typeof(udist).stringof);

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
            shuffle(proportions, rng);
        }
        watch.stop();
        writeln("\t", watch.peek.msecs, " ms. to generate ", l, " in-place random shuffles of ",
                proportions.length, " elements.");

        writeln();
    }
}
