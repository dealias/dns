TRI = $(HOME)/tri
ARCH = unix
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
vpath %.cc $(HOME)/fftw++
INCL = -I. -I$(TRI) -I$(HOME)/nw -I$(HOME)/fftw++

FILES = dns fftw++ convolution $(CORE) $(UTILS)
LIB += -lfftw3

include $(TRI)/config/Rules

