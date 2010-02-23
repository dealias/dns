all: conv cconv conv2 cconv2

#For g++
CC=g++
CFLAGS=-g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse -march=native

#For icpc
#CC=icpc
#CFLAGS=-O3 -ansi-alias -malign-double -fp-model fast=2 

LDFLAGS=-lfftw3

icpc: all
	CC=icpc

conv: conv.cc fftw++.cc
	$(CC)  $(CFLAGS) conv.cc fftw++.cc $(LDFLAGS) -o conv

cconv: cconv.cc fftw++.cc
	$(CC)  $(CFLAGS) cconv.cc fftw++.cc $(LDFLAGS) -o cconv

conv2: conv2.cc fftw++.cc
	$(CC)  $(CFLAGS) conv2.cc fftw++.cc $(LDFLAGS) -o conv2

cconv2: cconv2.cc fftw++.cc
	$(CC)  $(CFLAGS) cconv2.cc fftw++.cc $(LDFLAGS) -o cconv2

clean: 
	rm -rf conv cconv conv2 cconv2
