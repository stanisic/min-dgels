CLAPACK=src/base
BLAS_LIB=$(CLAPACK)/blas_LINUX.a
F2C_LIB=$(CLAPACK)/F2CLIBS/libf2c.a
ADDITIONAL=src/additional

all: clean install run

run:
	gcc $(ADDITIONAL)/*.c main.c -lm $(GSL_LIB) $(BLAS_LIB) $(F2C_LIB) -I$(ADDITIONAL) -omain
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
