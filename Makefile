TRI = $(HOME)/tri
ARCH = unix
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
INCL = -I. -I$(TRI) -I$(HOME)/nw -I-

FILES = dns Cartesian $(CORE) $(UTILS)
LIB += -lfftw3

include $(TRI)/config/Rules

