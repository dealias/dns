all: conv cconv conv2 cconv2

CFLAGS += -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse fftw++.cc -march=native
LDFLAGS += -lfftw3 

conv: conv.cc fftw++.cc
cconv: cconv.cc fftw++.cc

conv2: conv2.cc fftw++.cc
cconv2: cconv2.cc fftw++.cc

clean: 
	rm -rf conv cconv conv2 cconv2
