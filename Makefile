TRI = $(HOME)/tri
ARCH = unix
FFT = fft
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
INCL = -I. -I$(TRI) -I$(HOME)/nw -I-

TRIAD = dns Cartesian rfft $(FFT) $(CORE) $(UTILS)

DEPEND = $(TRIAD)

include $(TRI)/config/Rules

