CLAPACK=src/base
BLAS_LIB=$(CLAPACK)/BUILD/libminblas.a
F2C_LIB=$(CLAPACK)/BUILD/libminf2c.a
ADDITIONAL=src/additional

all: clean install run

run:
	gcc $(ADDITIONAL)/*.c main.c -lm $(BLAS_LIB) $(F2C_LIB) -I$(ADDITIONAL) -omain
	./main

install:
	cd $(CLAPACK) && $(MAKE) -j5 blaslib
	cd $(CLAPACK) && $(MAKE) -j5 f2clib

clean:
	cd $(CLAPACK) && $(MAKE) clean
	rm -rf main


archive: clean
	rm -rf min-dgels-src.tar.gz
	tar -zcvf min-dgels-src.tar.gz src/
