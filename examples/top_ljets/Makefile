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
KLFITTERINCLUDE = $(KLFITTERDIR)/library/include
KLFITTERINCLUDE_EXTRAS = $(KLFITTERDIR)/extras/include
KLFITTERLIBS    = -L$(KLFITTERDIR) -lKLFitter -lKLFitterExtras
ifneq ($(MACHINE), $(MAC))
	KLFITTERLIBSO = $(KLFITTERDIR)/libKLFitter.so $(KLFITTERDIR)/libKLFitterExtras.so
else
	KLFITTERLIBSO = $(KLFITTERDIR)/libKLFitter.dylib $(KLFITTERDIR)/libKLFitterExtras.dylib
endif

# Programs
CXX  = g++
RM   = rm -f
ECHO = echo
LN   = ln -s

# flags and libs
CXXFLAGS += $(ROOTCFLAGS) $(BATCFLAGS) -I$(KLFITTERINCLUDE) -I$(KLFITTERINCLUDE_EXTRAS)
CXXFLAGS += -Wall -Wno-deprecated -O2 -ggdb
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
	LIBLINK = libKLFitterExtras.dylib
endif

$(CXXOBJS) : $(CXXSRCS) $(KLFITTERLIBSO)
	cd $(KLFITTERDIR) && $(MAKE) extras
	$(CXX) $(CXXFLAGS) $(LIBS) $< -o $@
	@if [ "$(MACHINE)" == "$(MAC)" ]; then if [ ! -e $(LIBLINK) ]; then $(LN) $(KLFITTERLIBSO) $(LIBLINK) ; fi; fi

$(KLFITTERLIBSO):

clean :
	$(RM) $(GARBAGE)
