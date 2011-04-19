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
KLFITTERLIBS    = -L$(KLFITTERDIR) -lKLFitter -L$(KLFITTERDIR) -lKLFitterExtras
ifneq ($(MACHINE), $(MAC))
	KLFITTERLIBSO = $(KLFITTERDIR)/libKLFitter.so
	KLFITTERLIBSO_EXTRAS = $(KLFITTERDIR)/libKLFitterExtras.so
else
	KLFITTERLIBSO = $(KLFITTERDIR)/libKLFitter.dylib
	KLFITTERLIBSO_EXTRAS = $(KLFITTERDIR)/libKLFitterExtras.dylib
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
CXXSRCS   = runKLFitter.cxx
CXXOBJS   = runKLFitter.exe
ifeq ($(MACHINE), $(MAC))
	LIBLINK = libKLFitter.dylib
	LIBLINK_EXTRAS = libKLFitterExtras.dylib
endif
GARBAGE   = $(CXXOBJS) $(LIBLINK) $(LIBLINK_EXTRAS)

$(CXXOBJS) : $(CXXSRCS) $(KLFITTERLIBSO) $(KLFITTERLIBSO_EXTRAS)
	$(CXX) $(CXXFLAGS) $(LIBS) $< -o $@
	@if [ "$(MACHINE)" == "$(MAC)" ]; then if [ ! -e $(LIBLINK) ]; then $(LN) $(KLFITTERLIBSO) $(LIBLINK) ; fi; fi
	@if [ "$(MACHINE)" == "$(MAC)" ]; then if [ ! -e $(LIBLINK_EXTRAS) ]; then $(LN) $(KLFITTERLIBSO_EXTRAS) $(LIBLINK_EXTRAS) ; fi; fi

$(KLFITTERLIBSO) : $(KLFITTERLIBSO_EXTRAS)

$(KLFITTERLIBSO_EXTRAS): update-lib

update-lib:
	@cd $(KLFITTERDIR) && $(MAKE) extras

clean :
	$(RM) $(GARBAGE)
