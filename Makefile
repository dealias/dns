TRI = $(HOME)/tri
ARCH = unix
FFT = fft
POLL = poll

include $(TRI)/config/Common

vpath %.cc $(HOME)/nw
INCL = -I. -I$(TRI) -I$(HOME)/nw -I-

FILES = dns Cartesian rfft $(FFT) $(CORE) $(UTILS)
LIB += $(LIBFFT)

include $(TRI)/config/Rules

