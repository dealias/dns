IDIR =$(HOME)/fftw++/

DEFS=-g

# GNU compiler
ifeq ($(shell $(CXX) -v 2>&1 | tail -n 1 | head -c 3),gcc)
CXXFLAGS=-O3 -fopenmp -g -Wall -ansi -DNDEBUG -fomit-frame-pointer \
	-fstrict-aliasing -ffast-math -msse2 -mfpmath=sse -march=native
#For valgrind:
#CXXFLAGS=-fopenmp -g -Wall -ansi -fomit-frame-pointer \
#	-fstrict-aliasing -ffast-math -msse2 -mfpmath=sse
endif

#Intel compiler
ifeq ($(shell $(CXX) -v 2>&1 | head -c 4),icpc)
CXXFLAGS=-O3 -openmp -ansi-alias -malign-double -fp-model fast=2
endif

#IBM compiler
ifeq ($(shell $(CXX) -qversion 2>&1 | head -c 3),IBM)
CXXFLAGS=-O5 -P -qsmp -qalign -qarch -qtune -qcache -qipa -qarch=qp
endif

CXXFLAGS += $(DEFS) -I$(IDIR) $(CFLAGS)

ifneq ($(strip $(FFTW_INCLUDE_PATH)),)
CXXFLAGS+=-I$(FFTW_INCLUDE_PATH)
endif

LDFLAGS=
ifneq ($(strip $(FFTW_LIB_PATH)),)
LDFLAGS+=-L$(FFTW_LIB_PATH)
endif

#LDFLAGS=-lfftw3_threads -lfftw3 -lm
LDFLAGS+=-lfftw3_omp -lfftw3 -lm

MAKEDEPEND=$(CXXFLAGS) -O0 -M -DDEPEND

vpath %.cc $(IDIR)

FILES=dns
FFTW=fftw++
EXTRA=$(FFTW) convolution
ALL=$(FILES) $(EXTRA)

all: $(FILES)

%.o : %.cc %.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

dns: dns.o $(EXTRA:=.o)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

clean:  FORCE
	rm -rf $(ALL) $(ALL:=.o) $(ALL:=.d)

.SUFFIXES: .c .cc .o .d

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCL) -o $@ -c $<
.cc.d:
	@echo Creating $@; \
	rm -f $@; \
	${CXX} $(MAKEDEPEND) $(INCL) $< > $@.$$$$ 2>/dev/null && \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifeq (,$(findstring clean,${MAKECMDGOALS}))
-include $(ALL:=.d)
endif

FORCE:
