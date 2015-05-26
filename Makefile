TRI = $(HOME)/tri
ARCH = unix
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
vpath %.cc $(HOME)/fftw++
INCL += -I$(HOME)/nw -I$(HOME)/fftw++ -I$(HOME)/fftw/include

EXTRA = dnsbase fftw++ convolution explicit $(CORE) $(UTILS)
FILES = dns $(EXTRA)
OTHER = 

LIB += -L$(HOME)/fftw/lib -lfftw3 -lfftw3_omp

include $(TRI)/config/Rules

all: dns mdns

dns: 	dependencies
	+make -f $(TRI)/config/Compile FILES="dns $(EXTRA)" NAME=dns
