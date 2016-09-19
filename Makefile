SOURCE=src
include $(SOURCE)/min-make.inc

all: clean install run

run:
	gcc main.c $(SOURCE)/$(BUILD)/*.a -lm -I$(SOURCE)/$(ADDITIONAL) -omain
	./main

install:
	cd $(SOURCE) && $(MAKE) install

clean:
	cd $(SOURCE) && $(MAKE) clean
	rm -rf main

archive: clean
	rm -rf min-dgels-src.tar.gz
	tar -zcvf min-dgels-src.tar.gz src/
