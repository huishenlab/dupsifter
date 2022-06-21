CC = gcc
CFLAGS = -W -Wall -finline-functions -fPIC -std=gnu99 -Wno-unused-result -O3
CLIB = -lpthread -lz -lm -llzma -lbz2 -lcurl
CF_OPTIMIZE = 1

OS := $(shell uname)
ifeq ($(OS),  Darwin)
	CFLAGS += -Wno-unused-function
else
	CLIB += -lrt -ltinfo
endif

INCLUDE = include

########################
### different modes ####
########################

PROG = dupsifter dupsifter_khashl

.PHONY : setdebug debug build

build: exportcf $(PROG)

debug: CF_OPTIMIZE := 0
debug: CFLAGS += -pg
#debug: CFLAGS += -g -pg
debug: CFLAGS := $(filter-out -O3,$(CFLAGS))
debug: build

exportcf:
	$(eval export CF_OPTIMIZE)

#####################
##### libraries #####
#####################

LHTSLIB_DIR = htslib-1.15.1
LHTSLIB_INCLUDE = $(LHTSLIB_DIR)/htslib
LHTSLIB = $(LHTSLIB_DIR)/libhts.a
$(LHTSLIB) :
	make -C $(LHTSLIB_DIR) libhts.a

LKLIB_DIR = klib
LKLIB = $(LKLIB_DIR)/klib2.a
$(LKLIB) :
	make -C $(LKLIB_DIR) klib2.a

# Main program
LIBS = $(LKLIB) $(LHTSLIB)
dupsifter: $(LIBS) dupsifter.o
	$(CC) $(CFLAGS) dupsifter.o -o $@ -I$(LKLIB_DIR) -I$(LHTSLIB_INCLUDE) $(LIBS) $(CLIB)

dupsifter.o: dupsifter.c
	$(CC) -c $(CFLAGS) dupsifter.c -o $@ -I$(LKLIB_DIR) -I$(LHTSLIB_INCLUDE)

dupsifter_khashl: $(LIBS) dupsifter_khashl.o
	$(CC) $(CFLAGS) dupsifter_khashl.o -o $@ -I$(LKLIB_DIR) -I$(LHTSLIB_INCLUDE) $(LIBS) $(CLIB)

dupsifter_khashl.o: dupsifter_khashl.c
	$(CC) -c $(CFLAGS) dupsifter_khashl.c -o $@ -I$(LKLIB_DIR) -I$(LHTSLIB_INCLUDE)

# Clean
.PHONY: clean
clean:
	rm -f dupsifter dupsifter.o
	rm -f dupsifter_khashl dupsifter_khashl.o

purge: clean
	make -C $(LKLIB_DIR) purge
	make -C $(LHTSLIB_DIR) clean
