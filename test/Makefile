SOURCE=../src
include $(SOURCE)/min-make.inc

GSL_LIB=-lgsl -lgslcblas

all: clean compare_gsl

compare_gsl:
	gcc compare_gsl.c $(SOURCE)/$(BUILD)/*.a $(GSL_LIB) -lm -I$(SOURCE)/$(ADDITIONAL) -ocompare_gsl
	./compare_gsl

clean:
	rm -rf compare_gsl *~
