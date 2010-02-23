all: conv cconv conv2 cconv2

#For g++
CFLAGS += -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse -march=native

#For icpc
#CFLAGS += -O3 -ansi-alias -malign-double -fp-model fast=2 

LDFLAGS += -lfftw3 

conv: conv.cc fftw++.cc
cconv: cconv.cc fftw++.cc

conv2: conv2.cc fftw++.cc
cconv2: cconv2.cc fftw++.cc

clean: 
	rm -rf conv cconv conv2 cconv2
