TRI = $(HOME)/tri
ARCH = unix
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
vpath %.cc $(HOME)/fftw++
INCL = -I. -I$(TRI) -I$(HOME)/nw -I$(HOME)/fftw++

EXTRA = fftw++ convolution $(CORE) $(UTILS)
FILES = $(EXTRA)
LIB += -lfftw3

include $(TRI)/config/Rules

dns: 
	make -f $(TRI)/config/Compile FILES="dns $(FILES)" NAME=dns

mdns: 
	make -f $(TRI)/config/Compile FILES="mdns $(FILES)" NAME=mdns
