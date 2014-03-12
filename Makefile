DC = dmd
SRC = std/random2/*.d

all: unit

unit: $(SRC)
	$(DC) -main -unittest -debug -cov -of$@ $^

doc: $(SRC)
	$(DC) -o- -Ddhtml $^ ../dlang.org/std.ddoc

.PHONY: clean

clean:
	rm -rf unit *.o *.di html
