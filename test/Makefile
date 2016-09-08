CLAPACK=../src/base
BLAS_LIB=$(CLAPACK)/BUILD/libminblas.a
F2C_LIB=$(CLAPACK)/BUILD/libminf2c.a
ADDITIONAL=../src/additional

GSL_LIB=-lgsl -lgslcblas

all: clean compare_gsl

compare_gsl:
	gcc $(ADDITIONAL)/*.c compare_gsl.c -lm $(GSL_LIB) $(BLAS_LIB) $(F2C_LIB) -I$(ADDITIONAL) -ocompare_gsl
	./compare_gsl

clean:
	rm -rf compare_gsl
