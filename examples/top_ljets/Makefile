# check if OS is linux (Linux) or mac (Darwin)
MACHINE = $(shell uname -s)
LINUX   = Linux
MAC     = Darwin

# Root variables
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs) -lMinuit

# BAT
BATCFLAGS = -I$(BATINSTALLDIR)
BATLIBS   = -L$(BATINSTALLDIR)/lib -lBAT

# KLFitter
KLFITTERDIR     = ../..
KLFITTERINCLUDE = $(KLFITTERDIR)/include
KLFITTERLIBS    = -L$(KLFITTERDIR) -lKLFitter
ifneq ($(MACHINE), $(MAC))
	KLFITTERLIBSO = $(KLFITTERDIR)/libKLFitter.so
else
	KLFITTERLIBSO = $(KLFITTERDIR)/libKLFitter.dylib
endif

# Programs
CXX  = g++
RM   = rm -f
ECHO = echo
LN   = ln -s

# flags and libs
CXXFLAGS += $(ROOTCFLAGS) $(BATCFLAGS) -I$(KLFITTERINCLUDE)
CXXFLAGS += -Wall -Wno-deprecated -O2 
ifneq ($(MACHINE), $(MAC))
	CXXFLAGS += -fPIC
endif
LIBS      = $(ROOTLIBS) $(BATLIBS) $(KLFITTERLIBS)

# files
CXXSRCS   = runKLFitter.c
CXXOBJS   = runKLFitter.exe
GARBAGE   = $(CXXOBJS)
ifeq ($(MACHINE), $(MAC))
	LIBLINK = libKLFitter.dylib
endif

$(CXXOBJS) : $(CXXSRCS) $(KLFITTERLIBSO)
	$(CXX) $(CXXFLAGS) $(LIBS) $< -o $@
	@if [ "$(MACHINE)" == "$(MAC)" ]; then if [ ! -e $(LIBLINK) ]; then $(LN) $(KLFITTERLIBSO) $(LIBLINK) ; fi; fi

$(KLFITTERLIBSO):

clean :
	$(RM) $(GARBAGE)

print :
	$(ECHO) compiler   : $(CXX)
	$(ECHO) c++ srcs   : $(CXXSRCS)
	$(ECHO) c++ objs   : $(CXXOBJS)
	$(ECHO) c++ flags  : $(CXXFLAGS)
	$(ECHO) libs       : $(LIBS)
	$(ECHO) machine    : $(MACHINE)
	$(ECHO) linux      : $(LINUX)
	$(ECHO) mac        : $(MAC)
