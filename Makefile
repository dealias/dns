TRI = $(HOME)/tri
ARCH = unix
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
vpath %.cc $(HOME)/fftw++
INCL = -I. -I$(TRI) -I$(HOME)/nw -I$(HOME)/fftw++

EXTRA = fftw++ convolution $(CORE) $(UTILS)
FILES = dns $(EXTRA)
OTHER = mdns 

LIB += -lfftw3

include $(TRI)/config/Rules

all: dns mdns

dns: 	dependencies
	+make -f $(TRI)/config/Compile FILES="dns $(EXTRA)" NAME=dns

mdns: 	dependencies
	+make -f $(TRI)/config/Compile FILES="mdns $(EXTRA)" NAME=mdns
