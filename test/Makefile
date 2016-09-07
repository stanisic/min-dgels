CLAPACK=../src/base
BLAS_LIB=$(CLAPACK)/blas_LINUX.a
F2C_LIB=$(CLAPACK)/F2CLIBS/libf2c.a
ADDITIONAL=../src/additional

GSL_LIB=-lgsl -lgslcblas

all: clean compare_gsl

compare_gsl:
	gcc $(ADDITIONAL)/*.c compare_gsl.c -lm $(GSL_LIB) $(BLAS_LIB) $(F2C_LIB) -I$(ADDITIONAL) -ocompare_gsl
	./compare_gsl

clean:
	rm -rf compare_gsl
