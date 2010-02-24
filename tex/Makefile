all: conv cconv conv2 cconv2

#For g++
CC=g++
CFLAGS=-g -Wall -ansi -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse -march=native

#For icpc
a:=$(shell which icpc 2>&1 | tail -c5)
ifeq ($(a),icpc)
CC=icpc
CFLAGS=-O3 -ansi-alias -malign-double -fp-model fast=2 
endif

OPTS = $(CFLAGS)
MAKEDEPEND = $(OPTS) -O0 -M -DDEPEND
LDFLAGS=-lfftw3

FILES=conv cconv conv2 cconv2 fftw++

icpc: all
	CC=icpc

conv: conv.o fftw++.o
	$(CC)  $(CFLAGS) conv.o fftw++.o $(LDFLAGS) -o conv

cconv: cconv.o fftw++.o
	$(CC)  $(CFLAGS) cconv.o fftw++.o $(LDFLAGS) -o cconv

conv2: conv2.o fftw++.o
	$(CC)  $(CFLAGS) conv2.o fftw++.o $(LDFLAGS) -o conv2

cconv2: cconv2.o fftw++.o
	$(CC)  $(CFLAGS) cconv2.o fftw++.o $(LDFLAGS) -o cconv2

clean:  FORCE
	rm -rf $(FILES) $(FILES:=.o) $(FILES:=.d)

.SUFFIXES: .c .cc .o .d
.cc.o:
	$(CXX) $(OPTS) $(INCL) -o $@ -c $<
.cc.d:
	@echo Creating $@; \
	rm -f $@; \
	${CXX} $(MAKEDEPEND) $(INCL) $(MDOPTS) $< > $@.$$$$ 2>/dev/null && \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifeq (,$(findstring clean,${MAKECMDGOALS}))
-include $(FILES:=.d)
endif

FORCE:
