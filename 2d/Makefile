TRI = $(HOME)/triad
ARCH = unix
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
vpath %.cc $(HOME)/fftw++ $(HOME)/fftwpp

INCL += -I$(HOME)/nw -I$(HOME)/fftw++ -I$(HOME)/fftwpp -I$(HOME)/fftw/include

ifneq ($(strip $(FFTW_INCLUDE_PATH)),)
INCL += -I$(FFTW_INCLUDE_PATH)
endif

ifneq ($(strip $(FFTW_LIB_PATH)),)
LIB += -L$(FFTW_LIB_PATH)
endif

EXTRA = fftw++ convolve parallel $(CORE) $(UTILS)
FILES = dns $(EXTRA)
OTHER =

LIB += -lfftw3 -lfftw3_omp

include $(TRI)/config/Rules

all: dns mdns

dns: 	dependencies
	+make -f $(TRI)/config/Compile FILES="dns $(EXTRA)" NAME=dns
