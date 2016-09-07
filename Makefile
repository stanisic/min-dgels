CLAPACK=src/base
BLAS_LIB=$(CLAPACK)/blas_LINUX.a
F2C_LIB=$(CLAPACK)/F2CLIBS/libf2c.a
ADDITIONAL=src/additional

all: clean install run

run: #Write main

install:
	cd $(CLAPACK) && $(MAKE) -j5 blaslib
	cd $(CLAPACK) && $(MAKE) -j5 f2clib

clean:
	cd $(CLAPACK) && $(MAKE) clean
