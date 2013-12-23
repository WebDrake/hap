module std.random2.adaptor;

import std.random2.generator, std.random2.traits;

import std.range, std.traits;

// dice
/**
Rolls a random die with relative probabilities stored in $(D proportions).
Returns the index in $(D proportions) that was chosen.

Example:
----
auto x = dice(0.5, 0.5);   // x is 0 or 1 in equal proportions
auto y = dice(50, 50);     // y is 0 or 1 in equal proportions
auto z = dice(70, 20, 10); // z is 0 70% of the time, 1 20% of the time,
                           // and 2 10% of the time
----
*/
size_t dice(RandomGen, Num)(ref RandomGen rng, Num[] proportions...)
    if (isNumeric!Num && isForwardRange!RandomGen)
{
    return diceImpl(rng, proportions);
}

/// Ditto
size_t dice(RandomGen, Range)(ref RandomGen rng, Range proportions)
    if (isUniformRNG!RandomGen && isForwardRange!Range &&
        isNumeric!(ElementType!Range) && !isArray!Range)
{
    return diceImpl(rng, proportions);
}

/// Ditto
size_t dice(Range)(Range proportions)
    if (isForwardRange!Range && isNumeric!(ElementType!Range) && !isArray!Range)
{
    return diceImpl(rndGen, proportions);
}

/// Ditto
size_t dice(Num)(Num[] proportions...)
    if (isNumeric!Num)
{
    return diceImpl(rndGen, proportions);
}

private size_t diceImpl(RandomGen, Range)(ref RandomGen rng, Range proportions)
    if (isUniformRNG!RandomGen && isForwardRange!Range && isNumeric!(ElementType!Range))
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

// RandomCover
/**
 * Covers a given range $(D r) in a random manner, i.e. goes through each
 * element of $(D r) once and only once, but in a random order.  $(D r)
 * must be a random-access range with length.
 *
 * If no random number generator is passed to $(D randomCover), the
 * thread-global RNG rndGen will be used.
 *
 * Example:
 * ----
 * int[] a = [ 0, 1, 2, 3, 4, 5, 6, 7, 8 ];
 * foreach (e; randomCover(a))
 * {
 *     writeln(e);
 * }
 *
 * // using a specified random number generator
 * auto gen = new Random(unpredictableSeed);
 * foreach (e; randomCover(a, gen))
 * {
 *     writeln(e);
 * }
 * ----
 */
final class RandomCover(Range, RandomGen)
    if (isRandomAccessRange!Range && isUniformRNG!RandomGen)
{
  private:
    bool[] _chosen;
    size_t _current;
    size_t _alreadyChosen;
    Range _input;
    RandomGen _rng;

    // private constructor for use by .save
    this(Range input, RandomGen rng, in bool[] chosen, in size_t current, in size_t already)
    {
        _input = input;
        _rng = rng;
        _chosen.length = chosen.length;
        _chosen[] = chosen[];
        _current = current;
        assert(_chosen[_current]);
        _alreadyChosen = already;
        assert(_alreadyChosen > 0);
    }

  public:
    this(Range input, RandomGen rng)
    {
        _input = input;
        _rng = rng;
        _chosen.length = _input.length;
        _alreadyChosen = 0;
        popFront();
    }

    /// Range primitives.
    bool empty() @property @safe const nothrow pure
    {
        return _alreadyChosen > _input.length;
    }

    /// ditto
    auto ref front() @property @safe const nothrow pure
    in
    {
        assert(_alreadyChosen > 0);
    }
    body
    {
        return _input[_current];
    }

    /// ditto
    void popFront()
    {
        if (_alreadyChosen >= _input.length)
        {
            // No more elements
            ++_alreadyChosen; // means we're done
            return;
        }

        size_t k = _input.length - _alreadyChosen;
        size_t i;

        foreach (e; _input)
        {
            if (_chosen[i])
            {
                ++i;
                continue;
            }

            // Roll a dice with k faces
            import std.random2.distribution;
            auto chooseMe = uniform(0, k, _rng) == 0;
            assert(k > 1 || chooseMe);

            if (chooseMe)
            {
                _chosen[i] = true;
                _current = i;
                ++_alreadyChosen;
                return;
            }

            --k;
            ++i;
        }
    }

    /// ditto
    static if (hasLength!Range)
    {
        size_t length() @property @safe const nothrow pure
        in
        {
            assert(_alreadyChosen > 0);
        }
        body
        {
            return (1 + _input.length) - _alreadyChosen;
        }
    }

    /// ditto
    static if (isForwardRange!RandomGen)
    {
        typeof(this) save() @property
        {
            auto ret = new typeof(this)(_input.save, _rng.save, _chosen,
                                        _current, _alreadyChosen);
            return ret;
        }
    }

}

/// Ditto
auto randomCover(Range, RandomGen)(Range r, RandomGen rng)
    if (isRandomAccessRange!Range && isUniformRNG!RandomGen)
{
    return new RandomCover!(Range, RandomGen)(r, rng);
}

/// Ditto
auto randomCover(Range)(Range r)
    if (isRandomAccessRange!Range)
{
    return new RandomCover!(Range, Random)(r, rndGen);
}

unittest
{
    import std.algorithm, std.string, std.typetuple;

    int[] a = [ 0, 1, 2, 3, 4, 5, 6, 7, 8 ];

    foreach (RandomGen; TypeTuple!(void, UniformRNGTypes))
    {
        static if (is(RandomGen == void))
        {
            auto rc = randomCover(a);
            static assert(isInputRange!(typeof(rc)));
            static assert((isForwardRange!Random && isForwardRange!(typeof(rc))) ||
                          (!isForwardRange!Random && !isForwardRange!(typeof(rc))));
        }
        else
        {
            auto rng = new RandomGen(unpredictableSeed);
            auto rc = randomCover(a, rng);
            static assert(isInputRange!(typeof(rc)));
            static assert((isForwardRange!RandomGen && isForwardRange!(typeof(rc))) ||
                          (!isForwardRange!RandomGen && !isForwardRange!(typeof(rc))));
        }

        auto rc2 = rc.save;

        int[] b = new int[9];
        uint i;
        foreach (e, e2; lockstep(rc, rc2))
        {
            b[i++] = e;
            assert(e == e2);
        }
        sort(b);
        assert(a == b, format("%s", b));
    }
}

// RandomSample
/**
 * Selects a random subsample out of $(D r), containing exactly $(D n)
 * elements. The order of elements is the same as in the original
 * range. The total length of $(D r) must be known. If $(D total) is
 * passed in, the total number of sample is considered to be $(D
 * total). Otherwise, $(D RandomSample) uses $(D r.length).
 *
 * $(D RandomSample) implements Jeffrey Scott Vitter's Algorithm D
 * (see Vitter $(WEB dx.doi.org/10.1145/358105.893, 1984), $(WEB
 * dx.doi.org/10.1145/23002.23003, 1987)), which selects a sample
 * of size $(D n) in O(n) steps and requiring O(n) random variates,
 * regardless of the size of the data being sampled.  The exception
 * to this is if traversing k elements on the input range is itself
 * an O(k) operation (e.g. when sampling lines from an input file),
 * in which case the sampling calculation will inevitably be of
 * O(total).
 *
 * RandomSample will throw an error if $(D total) is verifiably less
 * than the total number of elements available in the input, or if
 * $(D n > total).
 *
 * If no random number generator is passed to $(D randomSample), the
 * thread-local default RNG rndGen will be used to generate the
 * sample.
 *
 * Example:
 * ----
 * int[] arr = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ];
 * // Print 5 random elements picked from arr
 * foreach (e; randomSample(arr, 5))
 * {
 *     writeln(e);
 * }
 *
 * // Print 5 random elements picked from arr,
 * // using a specified random number generator
 * auto gen = new Random(unpredictableSeed);
 * foreach (e; randomSample(arr, 5, gen))
 * {
 *     writeln(e);
 * }
 * ----
 */
final class RandomSample(Range, RandomGen)
    if (isInputRange!Range && isUniformRNG!RandomGen)
{
  private:
    import std.random2.distribution;
    enum ushort _alphaInverse = 13; // Vitter's recommended value.
    enum Skip { None, A, D };
    size_t _available, _toSelect, _index;
    double _Vprime;
    Range _input;
    Skip _skip = Skip.None;
    RandomGen _rng;

    // private constructor, for use with .save
    this(Range input, RandomGen rng)
    {
        _input = input;
        _rng = rng;
    }

    // Randomly reset the value of _Vprime.
    double newVprime(size_t remaining)
    {
        return uniform!"()"(0.0, 1.0, _rng) ^^ (1.0 / remaining);
    }

    void prime()
    {
        if (empty)
        {
            return;
        }
        assert(_available && _available >= _toSelect);
        immutable size_t s = skip();
        assert(s + _toSelect <= _available);
        static if (hasLength!Range)
        {
            assert(s + _toSelect <= _input.length);
        }
        assert(!_input.empty);
        _input.popFrontExactly(s);
        _index += s;
        _available -= s;
        assert(_available > 0);
    }

    size_t skip()
    in
    {
        assert(_skip != Skip.None);
    }
    body
    {
        // Step D1: if the number of points still to select is greater
        // than a certain proportion of the remaining data points, i.e.
        // if n >= alpha * N where alpha = 1/13, we carry out the
        // sampling with Algorithm A.
        if (_skip == Skip.A)
        {
            return skipA();
        }
        else if ((_alphaInverse * _toSelect) > _available)
        {
            // We shouldn't get here unless the current selected
            // algorithm is D.
            assert(_skip == Skip.D);
            _skip = Skip.A;
            return skipA();
        }
        else
        {
            assert(_skip == Skip.D);
            return skipD();
        }
    }

    /* Vitter's Algorithm A, used when the ratio of needed sample values
     * to remaining data values is sufficiently large.
     */
    size_t skipA()
    {
        size_t s;
        double v, quot, top;

        if (_toSelect==1)
        {
            s = uniform(0, _available, _rng);
        }
        else
        {
            top = _available - _toSelect;
            quot = top / _available;
            v = uniform!"()"(0.0, 1.0, _rng);

            while (quot > v)
            {
                ++s;
                quot *= (top - s) / (_available - s);
            }
        }

        return s;
    }

    /* Vitter's Algorithm D.  For an extensive description of the algorithm
     * and its rationale, see:
     *
     *    * Vitter, J.S. (1984), "Faster methods for random sampling",
     *      Commun. ACM 27(7): 703--718
     *
     *    * Vitter, J.S. (1987) "An efficient algorithm for sequential
     *      random sampling", ACM Trans. Math. Softw. 13(1): 58-67.
     *
     * Variable names are chosen to match those in Vitter's paper.
     */
    size_t skipD()
    {
        import std.math;

        // Confirm that the check in Step D1 is valid and we
        // haven't been sent here by mistake
        assert((_alphaInverse * _toSelect) <= _available);

        // Now it's safe to use the standard Algorithm D mechanism.
        if (_toSelect > 1)
        {
            size_t s;
            size_t qu1 = 1 + _available - _toSelect;
            double x, y1;

            assert(!_Vprime.isNaN);

            while (true)
            {
                // Step D2: set values of x and u.
                for (x = _available * (1-_Vprime), s = cast(size_t) trunc(x);
                     s >= qu1;
                     x = _available * (1-_Vprime), s = cast(size_t) trunc(x))
                {
                    _Vprime = newVprime(_toSelect);
                }

                double u = uniform!"()"(0.0, 1.0, _rng);

                y1 = (u * (cast(double) _available) / qu1) ^^ (1.0/(_toSelect - 1));

                _Vprime = y1 * ((-x/_available)+1.0) * ( qu1/( (cast(double) qu1) - s ) );

                // Step D3: if _Vprime <= 1.0 our work is done and we return S.
                // Otherwise ...
                if (_Vprime > 1.0)
                {
                    size_t top = _available - 1, limit;
                    double y2 = 1.0, bottom;

                    if (_toSelect > (s+1))
                    {
                        bottom = _available - _toSelect;
                        limit = _available - s;
                    }
                    else
                    {
                        bottom = _available - (s+1);
                        limit = qu1;
                    }

                    foreach (size_t t; limit .. _available)
                    {
                        y2 *= top/bottom;
                        top--;
                        bottom--;
                    }

                    // Step D4: decide whether or not to accept the current value of S.
                    if (_available/(_available-x) < y1 * (y2 ^^ (1.0/(_toSelect-1))))
                    {
                        // If it's not acceptable, we generate a new value of _Vprime
                        // and go back to the start of the for (;;) loop.
                        _Vprime = newVprime(_toSelect);
                    }
                    else
                    {
                        // If it's acceptable we generate a new value of _Vprime
                        // based on the remaining number of sample points needed,
                        // and return S.
                        _Vprime = newVprime(_toSelect-1);
                        return s;
                    }
                }
                else
                {
                    // Return if condition D3 satisfied.
                    return s;
                }
            }
        }
        else
        {
            // If only one sample point remains to be taken ...
            return cast(size_t) trunc(_available * _Vprime);
        }
    }

  public:
    static if (hasLength!Range)
    {
        this(Range input, size_t howMany, RandomGen rng)
        {
            this(input, howMany, input.length, rng);
        }
    }

    this(Range input, size_t howMany, size_t total, RandomGen rng)
    {
        import std.exception, std.string : format;
        _input = input;
        _rng = rng;
        _available = total;
        _toSelect = howMany;
        enforce(_toSelect <= _available,
                format("RandomSample: cannot sample %s items when only %s are available",
                       _toSelect, _available));
        static if (hasLength!Range)
        {
            enforce(_available <= _input.length,
                    format("RandomSample: specified %s items as available when input contains only %s",
                           _available, _input.length));
        }

        /* We can save ourselves a random variate by checking right
         * at the beginning if we should use Algorithm A.
         */
        if ((_alphaInverse * _toSelect) > _available)
        {
            _skip = Skip.A;
        }
        else
        {
            _skip = Skip.D;
            _Vprime = newVprime(_toSelect);
        }
        prime();
    }


    /// Range primitives.
    bool empty() @property @safe const nothrow pure
    {
        return _toSelect == 0;
    }

    /// ditto
    auto ref front() @property
    in
    {
        /* Check if sample has been initialized and that
         * there are still sample points left to generate
         */
        assert(_skip != Skip.None);
        assert(!empty);
    }
    body
    {
        return _input.front;
    }

    /// ditto
    void popFront()
    in
    {
        // Check that sample has been initialized
        assert(_skip != Skip.None);
    }
    body
    {
        _input.popFront();
        --_available;
        --_toSelect;
        ++_index;
        prime();
    }

    /// ditto
    static if (isForwardRange!Range && isForwardRange!RandomGen)
    {
        typeof(this) save() @property
        {
            auto ret = new typeof(this)(_input.save, _rng.save);
            ret._available = this._available;
            ret._toSelect = this._toSelect;
            ret._Vprime = this._Vprime;
            ret._index = this._index;
            ret._skip = this._skip;
            return ret;
        }
    }

    /// ditto
    size_t length() @property @safe const nothrow pure
    {
        return _toSelect;
    }

    /// Returns the index of the visited record.
    size_t index()
    in
    {
        assert(_skip != Skip.None);
    }
    body
    {
        return _index;
    }

}

/// Ditto
auto randomSample(Range)(Range r, size_t n, size_t total)
    if (isInputRange!Range)
{
    return new RandomSample!(Range, Random)(r, n, total, rndGen);
}

/// Ditto
auto randomSample(Range)(Range r, size_t n)
    if (isInputRange!Range && hasLength!Range)
{
    return new RandomSample!(Range, Random)(r, n, r.length, rndGen);
}

/// Ditto
auto randomSample(Range, RandomGen)(Range r, size_t n, size_t total, RandomGen rng)
    if (isInputRange!Range && isUniformRNG!RandomGen)
{
    return new RandomSample!(Range, RandomGen)(r, n, total, rng);
}

/// Ditto
auto randomSample(Range, RandomGen)(Range r, size_t n, RandomGen rng)
    if (isInputRange!Range && hasLength!Range && isUniformRNG!RandomGen)
{
    return new RandomSample!(Range, RandomGen)(r, n, r.length, rng);
}

unittest
{
    // For test purposes, an infinite input range
    struct TestInputRange
    {
        private auto r = recurrence!"a[n-1] + 1"(0);
        bool empty() @property const pure nothrow { return r.empty; }
        auto front() @property pure nothrow { return r.front; }
        void popFront() pure nothrow { r.popFront(); }
    }
    static assert(isInputRange!TestInputRange);
    static assert(!isForwardRange!TestInputRange);

    int[] a = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ];

    foreach (UniformRNG; UniformRNGTypes)
    {
        auto rng = new UniformRNG(unpredictableSeed);
        /* First test the most general case: randomSample of input range, with and
         * without a specified random number generator.
         */
        static assert(isInputRange!(typeof(randomSample(TestInputRange(), 5, 10))));
        static assert(isInputRange!(typeof(randomSample(TestInputRange(), 5, 10, rng))));
        static assert(!isForwardRange!(typeof(randomSample(TestInputRange(), 5, 10))));
        static assert(!isForwardRange!(typeof(randomSample(TestInputRange(), 5, 10, rng))));
        // test case with range initialized by direct call to struct
        {
            auto sample =
                new RandomSample!(TestInputRange, UniformRNG)
                                 (TestInputRange(), 5, 10, rng);
            static assert(isInputRange!(typeof(sample)));
            static assert(!isForwardRange!(typeof(sample)));
            assert(sample.length == 5);
            assert(sample._available == 10 - sample.index);
        }

        /* Now test the case of an input range with length.  We ignore the cases
         * already covered by the previous tests.
         */
        static assert(isInputRange!(typeof(randomSample(TestInputRange().takeExactly(10), 5))));
        static assert(isInputRange!(typeof(randomSample(TestInputRange().takeExactly(10), 5, rng))));
        static assert(!isForwardRange!(typeof(randomSample(TestInputRange().takeExactly(10), 5))));
        static assert(!isForwardRange!(typeof(randomSample(TestInputRange().takeExactly(10), 5, rng))));
        // test case with range initialized by direct call to struct
        {
            auto sample =
                new RandomSample!(typeof(TestInputRange().takeExactly(10)), UniformRNG)
                                 (TestInputRange().takeExactly(10), 5, 10, rng);
            static assert(isInputRange!(typeof(sample)));
            static assert(!isForwardRange!(typeof(sample)));
            assert(sample.length == 5);
            assert(sample._available == 10 - sample.index);
        }

        // Now test the case of providing a forward range as input.
        static assert(isForwardRange!(typeof(randomSample(a, 5))));
        static if (isForwardRange!UniformRNG)
        {
            static assert(isForwardRange!(typeof(randomSample(a, 5, rng))));
            // ... and test with range initialized directly
            {
                auto sample = new RandomSample!(int[], UniformRNG)(a, 5, rng);
                static assert(isForwardRange!(typeof(sample)));
                assert(sample.length == 5);
                assert(sample._available == a.length - sample.index);
            }
        }
        else
        {
            static assert(isInputRange!(typeof(randomSample(a, 5, rng))));
            static assert(!isForwardRange!(typeof(randomSample(a, 5, rng))));
            // ... and test with range initialized directly
            {
                auto sample = new RandomSample!(int[], UniformRNG)(a, 5, rng);
                static assert(isInputRange!(typeof(sample)));
                static assert(!isForwardRange!(typeof(sample)));
                assert(sample.length == 5);
                assert(sample._available == a.length - sample.index);
            }
        }

        /* Check that randomSample will throw an error if we claim more
         * items are available than there actually are, or if we try to
         * sample more items than are available. */
        import std.exception;
        assert(collectExceptionMsg(randomSample(a, 5, 15)) == "RandomSample: specified 15 items as available when input contains only 10");
        assert(collectExceptionMsg(randomSample(a, 15)) == "RandomSample: cannot sample 15 items when only 10 are available");
        assert(collectExceptionMsg(randomSample(a, 9, 8)) == "RandomSample: cannot sample 9 items when only 8 are available");
        assert(collectExceptionMsg(randomSample(TestInputRange(), 12, 11)) == "RandomSample: cannot sample 12 items when only 11 are available");

        /* Check that sampling algorithm never accidentally overruns the end of
         * the input range.  If input is an InputRange without .length, this
         * relies on the user specifying the total number of available items
         * correctly.
         */
        {
            uint i = 0;
            foreach (e; randomSample(a, a.length))
            {
                assert(e == i);
                ++i;
            }
            assert(i == a.length);

            i = 0;
            foreach (e; randomSample(TestInputRange(), 17, 17))
            {
                assert(e == i);
                ++i;
            }
            assert(i == 17);
        }


        // Check length properties of random samples.
        assert(randomSample(a, 5).length == 5);
        assert(randomSample(a, 5, 10).length == 5);
        assert(randomSample(a, 5, rng).length == 5);
        assert(randomSample(a, 5, 10, rng).length == 5);
        assert(randomSample(TestInputRange(), 5, 10).length == 5);
        assert(randomSample(TestInputRange(), 5, 10, rng).length == 5);

        // ... and emptiness!
        assert(randomSample(a, 0).empty);
        assert(randomSample(a, 0, 5).empty);
        assert(randomSample(a, 0, rng).empty);
        assert(randomSample(a, 0, 5, rng).empty);
        assert(randomSample(TestInputRange(), 0, 10).empty);
        assert(randomSample(TestInputRange(), 0, 10, rng).empty);

        /* Test that the (lazy) evaluation of random samples works correctly.
         *
         * We cover 2 different cases: a sample where the ratio of sample points
         * to total points is greater than the threshold for using Algorithm, and
         * one where the ratio is small enough (< 1/13) for Algorithm D to be used.
         *
         * For each, we also cover the case with and without a specified RNG.
         */
        {
            uint i = 0;

            // Small sample/source ratio, no specified RNG.
            foreach (e; randomSample(randomCover(a), 5))
            {
                ++i;
            }
            assert(i == 5);

            // Small sample/source ratio, specified RNG.
            i = 0;
            foreach (e; randomSample(randomCover(a), 5, rng))
            {
                ++i;
            }
            assert(i == 5);

            // Large sample/source ratio, no specified RNG.
            i = 0;
            foreach (e; randomSample(TestInputRange(), 123, 123_456))
            {
                ++i;
            }
            assert(i == 123);

            // Large sample/source ratio, specified RNG.
            i = 0;
            foreach (e; randomSample(TestInputRange(), 123, 123_456, rng))
            {
                ++i;
            }
            assert(i == 123);

            /* Sample/source ratio large enough to start with Algorithm D,
             * small enough to switch to Algorithm A.
             */
            i = 0;
            foreach (e; randomSample(TestInputRange(), 10, 131))
            {
                ++i;
            }
            assert(i == 10);
        }

        // Test that the .index property works correctly
        {
            auto sample1 = randomSample(TestInputRange(), 654, 654_321);
            for (; !sample1.empty; sample1.popFront())
            {
                assert(sample1.front == sample1.index);
            }

            auto sample2 = randomSample(TestInputRange(), 654, 654_321, rng);
            for (; !sample2.empty; sample2.popFront())
            {
                assert(sample2.front == sample2.index);
            }

            /* Check that it also works if .index is called before .front.
             * See: http://d.puremagic.com/issues/show_bug.cgi?id=10322
             */
            auto sample3 = randomSample(TestInputRange(), 654, 654_321);
            for (; !sample3.empty; sample3.popFront())
            {
                assert(sample3.index == sample3.front);
            }

            auto sample4 = randomSample(TestInputRange(), 654, 654_321, rng);
            for (; !sample4.empty; sample4.popFront())
            {
                assert(sample4.index == sample4.front);
            }
        }

        /* Test behaviour if .popFront() is called before sample is read.
         * This is a rough-and-ready check that the statistical properties
         * are in the ballpark -- not a proper validation of statistical
         * quality!  This incidentally also checks for reference-type
         * initialization bugs, as the foreach() loop will operate on a
         * copy of the popFronted (and hence initialized) sample.
         */
        {
            size_t count1, count99;
            foreach(_; 0 .. 100_000)
            {
                auto sample = randomSample(iota(100), 5);
                sample.popFront();
                foreach(s; sample)
                {
                    /* This is a sequential sampling process: 0 can only be
                     * the first sample point, so _can't_ be in the remainder
                     * of the sample after .popFront() is called.
                     */
                    assert(s != 0);

                    if (s == 1)
                    {
                        ++count1;
                    }
                    else if (s == 99)
                    {
                        ++count99;
                    }
                }
            }
            /* We already established that because this is a sequential sampling
             * process, 0 can only be the _first_ sample point and not one of the
             * sample points remaining after .popFront() is first called.  By the
             * same argument, 1 can only be in the remainder if it's the 2nd point
             * of the whole sample, and hence if 0 was the first; probability of 0
             * being first and 1 second is 5/100 * 4/99 (thank you, Algorithm S:-)
             * and so the mean count of 1 should be about 202.

             * Conversely, 99 can only be the _last_ sample point to be picked,
             * so its probability of inclusion should be independent of the initial
             * .popFront() and it should occur with frequency 5/100, hence its count
             * should be about 5000.
             *
             * Ideally we would perform far more trials and have a much higher
             * tolerance for this unittest, but the time required to do so is
             * out of all proportion to all other Phobos unittests, so we have
             * to put up with quite high variance in the outcome.  The current
             * test should at least highlight any extreme biases and could be
             * backed up by a hardcore external random test suite.
             */
            import std.string : format;
            assert(count1 < 300, format("1: %s > 300.", count1));
            assert(4_700 < count99, format("99: %s <= 4700.", count99));
            assert(count99 < 5_300, format("99: %s >= 5300.", count99));
        }

        /* Test that the save property works where input is a forward range,
         * and RandomSample is using a (forward range) random number generator
         * that is not rndGen.
         */
        static if (isForwardRange!UniformRNG)
        {
            auto sample1 = randomSample(a, 5, rng);
            auto sample2 = sample1.save;
            assert(sample1.array() == sample2.array());
        }

        // Bugzilla 8314
        {
            UniformRNG gen = new UniformRNG;

            auto sample(RandomGen)(uint seed, RandomGen gen)
                if (isUniformRNG!RandomGen)
            {
                gen.seed(seed);
                return randomSample(a, 1, gen).front;
            }

            // Start from 1 because not all RNGs accept 0 as seed.
            immutable fst = sample(1, gen);
            uint n = 1;
            while (sample!UniformRNG(++n, gen) == fst && n < n.max) {}
            assert(n < n.max);
        }
    }
}

/**
 * Shuffles elements of $(D r) using $(D gen) as a shuffler. $(D r) must be
 * a random-access range with length.  If no RNG is specified, $(D rndGen)
 * will be used.
 */
void randomShuffle(Range, RandomGen)(Range r, ref RandomGen gen)
    if(isRandomAccessRange!Range && isUniformRNG!RandomGen)
{
    return partialShuffle!(Range, RandomGen)(r, r.length, gen);
}

/// ditto
void randomShuffle(Range)(Range r)
    if(isRandomAccessRange!Range)
{
    return randomShuffle(r, rndGen);
}

unittest
{
    foreach(RandomGen; UniformRNGTypes)
    {
        // Also tests partialShuffle indirectly.
        auto a = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        auto b = a.dup;
        auto gen = new RandomGen(unpredictableSeed);
        randomShuffle(a, gen);
        assert(a.sort == b);
        randomShuffle(a);
        assert(a.sort == b);
    }
}

/**
 * Partially shuffles the elements of $(D r) such that upon returning $(D r[0..n])
 * is a random subset of $(D r) and is randomly ordered.  $(D r[n..r.length])
 * will contain the elements not in $(D r[0..n]).  These will be in an undefined
 * order, but will not be random in the sense that their order after
 * $(D partialShuffle) returns will not be independent of their order before
 * $(D partialShuffle) was called.
 *
 * $(D r) must be a random-access range with length.  $(D n) must be less than
 * or equal to $(D r.length).  If no RNG is specified, $(D rndGen) will be used.
 */
void partialShuffle(Range, RandomGen)(Range r, in size_t n, ref RandomGen gen)
    if(isRandomAccessRange!Range && isUniformRNG!RandomGen)
{
    import std.algorithm, std.exception, std.random2.distribution;
    enforce(n <= r.length, "n must be <= r.length for partialShuffle.");
    foreach (i; 0 .. n)
    {
        swapAt(r, i, i + uniform(0, n - i, gen));
    }
}

/// ditto
void partialShuffle(Range)(Range r, in size_t n)
    if(isRandomAccessRange!Range)
{
    return partialShuffle(r, n, rndGen);
}

unittest
{
    foreach(RandomGen; UniformRNGTypes)
    {
        auto a = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        auto b = a.dup;
        auto gen = new RandomGen(unpredictableSeed);
        partialShuffle(a, 5, gen);
        assert(a[5 .. $] == b[5 .. $]);
        assert(a[0 .. 5].sort == b[0 .. 5]);
        partialShuffle(a, 6);
        assert(a[6 .. $] == b[6 .. $]);
        assert(a[0 .. 6].sort == b[0 .. 6]);
    }
}
