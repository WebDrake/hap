DC = dmd

SRC = hap/random/adaptor.d      \
      hap/random/distribution.d \
      hap/random/generator.d    \
      hap/random/package.d      \
      hap/random/traits.d

EXPERIMENTAL = hap/random/device.d

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
	$(DC) -o- -Ddhtml $^ doc.ddoc

.PHONY: clean

clean:
	rm -rf benchmarknew benchmarkold unit unit-xper *.o *.di html
