TRI = $(HOME)/triad
ARCH = unix
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
vpath %.cc $(HOME)/fftw++
INCL += -I$(HOME)/nw -I$(HOME)/fftw++ -I$(HOME)/fftw/include

vpath %.cc $(HOME)/fftw++/tests # TEMPORARY until compact=false upgrade
INCL += -I$(HOME)/fftw++/tests # TEMPORARY until compact=false upgrade

ifneq ($(strip $(FFTW_INCLUDE_PATH)),)
INCL += -I$(FFTW_INCLUDE_PATH)
endif

ifneq ($(strip $(FFTW_LIB_PATH)),)
LIB += -L$(FFTW_LIB_PATH)
endif

EXTRA = fftw++ convolution explicit $(CORE) $(UTILS)
FILES = dns $(EXTRA)
OTHER = 

LIB += -L$(HOME)/fftw/lib -lfftw3 -lfftw3_omp

include $(TRI)/config/Rules

all: dns mdns

dns: 	dependencies
	+make -f $(TRI)/config/Compile FILES="dns $(EXTRA)" NAME=dns
