CC=g++
#CFLAGS=-g -Wall -ansi -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse -march=native
CFLAGS=-g -Wall -ansi -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -msse2 -mfpmath=sse -march=native

a:=$(shell which icpc 2>&1 | tail -c5)
ifeq ($(a),icpc)
CC=icpc
CFLAGS=-O3 -ansi-alias -malign-double -fp-model fast=2 
endif

MAKEDEPEND=$(CFLAGS) -O0 -M -DDEPEND
LDFLAGS=-lfftw3

FILES=conv cconv conv2 cconv2 conv3 cconv3 biconv biconv2
EXTRA=fftw++ convolution
ALL=$(FILES) $(EXTRA)

all: $(FILES)

icpc: all
	CC=icpc

conv: conv.o $(EXTRA:=.o)
	$(CC) $(CFLAGS) conv.o $(EXTRA:=.o) $(LDFLAGS) -o conv

cconv: cconv.o $(EXTRA:=.o)
	$(CC) $(CFLAGS) cconv.o $(EXTRA:=.o) $(LDFLAGS) -o cconv

conv2: conv2.o $(EXTRA:=.o)
	$(CC) $(CFLAGS) conv2.o $(EXTRA:=.o) $(LDFLAGS) -o conv2

cconv2: cconv2.o $(EXTRA:=.o)
	$(CC) $(CFLAGS) cconv2.o $(EXTRA:=.o) $(LDFLAGS) -o cconv2

conv3: conv3.o $(EXTRA:=.o)
	$(CC) $(CFLAGS) conv3.o $(EXTRA:=.o) $(LDFLAGS) -o conv3

cconv3: cconv3.o $(EXTRA:=.o)
	$(CC) $(CFLAGS) cconv3.o $(EXTRA:=.o) $(LDFLAGS) -o cconv3

biconv: biconv.o $(EXTRA:=.o)
	$(CC) $(CFLAGS) biconv.o $(EXTRA:=.o) $(LDFLAGS) -o biconv

biconv2: biconv2.o $(EXTRA:=.o)
	$(CC) $(CFLAGS) biconv2.o $(EXTRA:=.o) $(LDFLAGS) -o biconv2

clean:  FORCE
	rm -rf $(ALL) $(ALL:=.o) $(ALL:=.d)

.SUFFIXES: .c .cc .o .d
.cc.o:
	$(CXX) $(CFLAGS) $(INCL) -o $@ -c $<
.cc.d:
	@echo Creating $@; \
	rm -f $@; \
	${CXX} $(MAKEDEPEND) $(INCL) $< > $@.$$$$ 2>/dev/null && \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifeq (,$(findstring clean,${MAKECMDGOALS}))
-include $(ALL:=.d)
endif

FORCE:
