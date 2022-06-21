CC     = gcc
AR     = ar
CFLAGS = -g -Wall

ifeq (1, $(CF_OPTIMIZE))
	CFLAGS += -O2
	CFLAGS := $(filter-out -g,$(CFLAGS))
endif

# leave out for conflict with htslib when compiled by OSX clang
BGZF_OBJ = \
	bgzf.o

KLIB_OBJS = \
	kexpr.o \
	khmm.o \
	kmath.o \
	knetfile.o \
	knhx.o \
	kopen.o \
	ksa.o \
	kson.o \
	kstring.o \
	ksw.o \
	kthread.o \
	kurl.o

build: klib.a

.c.o :
	$(CC) -c $(CFLAGS) $< -o $@

klib.a: $(KLIB_OBJS) $(BGZF_OBJ)
	@-rm -f $@
	$(AR) -csr $@ $^

# without BGZF objects
klib2.a: $(KLIB_OBJS)
	@-rm -f $@
	$(AR) -csr $@ $^

clean:
	rm -f *.o **/*.o

purge: clean
	rm -f *.a **/*.a

cleanhist:
	rm -rf .git .gitignore
