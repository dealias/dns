TRI = $(HOME)/tri
ARCH = unix
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
INCL = -I. -I$(TRI) -I$(HOME)/nw

FILES = dns fftw++ convolution $(CORE) $(UTILS)
LIB += -lfftw3

include $(TRI)/config/Rules

