DC = dmd

SRC = source/hap/random/adaptor.d      \
      source/hap/random/distribution.d \
      source/hap/random/generator.d    \
      source/hap/random/package.d      \
      source/hap/random/traits.d

EXPERIMENTAL = source/hap/random/device.d

all: benchmark unit

benchmark: benchmarknew benchmarkold

benchmarknew: benchmarknew.d $(SRC)
	$(DC) -I./source -O -inline -release -of$@ $^

benchmarkold: benchmarkold.d $(SRC)
	$(DC) -I./source -O -inline -release -of$@ $^

unit: $(SRC)
	$(DC) -I./source -main -unittest -debug -cov -of$@ $^
	./unit

unit-xper: $(EXPERIMENTAL) $(SRC)
	$(DC) -I./source -main -unittest -debug -cov -of$@ $^
	./unit-xper

unittest: unit unit-xper

doc: $(SRC) $(EXPERIMENTAL)
	$(DC) -I./source -o- -Ddhtml $^ doc.ddoc template.ddoc

.PHONY: clean

clean:
	rm -rf benchmarknew benchmarkold unit unit-xper *.o *.di html
