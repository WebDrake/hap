DC = dmd
SRC = std/random2/*.d

all: benchmark unit

benchmark: benchmarknew

benchmarknew: benchmarknew.d $(SRC)
	$(DC) -O -inline -release -of$@ $^

unit: $(SRC)
	$(DC) -main -unittest -debug -cov -of$@ $^

doc: $(SRC)
	$(DC) -o- -Ddhtml $^ ../dlang.org/std.ddoc

.PHONY: clean

clean:
	rm -rf benchmarknew unit *.o *.di html
