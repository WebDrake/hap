DC = dmd

SRC = std/random2/adaptor.d      \
      std/random2/distribution.d \
      std/random2/generator.d    \
      std/random2/package.d      \
      std/random2/traits.d

EXPERIMENTAL = std/random2/device.d

all: benchmark unit

benchmark: benchmarknew benchmarkold

benchmarknew: benchmarknew.d $(SRC)
	$(DC) -O -inline -release -of$@ $^

benchmarkold: benchmarkold.d $(SRC)
	$(DC) -O -inline -release -of$@ $^

unit: $(SRC)
	$(DC) -main -unittest -debug -cov -of$@ $^

unit-xper: $(EXPERIMENTAL) $(SRC)
	$(DC) -main -unittest -debug -cov -of$@ $^

unittest: unit unit-xper
	./unit
	./unit-xper

doc: $(SRC) $(EXPERIMENTAL)
	$(DC) -o- -Ddhtml $^ ../dlang.org/std.ddoc

.PHONY: clean

clean:
	rm -rf benchmarknew benchmarkold unit unit-xper *.o *.di html
