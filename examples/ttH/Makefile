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
KLFITTERDIR = ../..
KLFITTERDIR_EXTRAS = $(KLFITTERDIR)/KLFitterExtras
KLFITTERLIBS = -L$(KLFITTERDIR)/lib -lKLFitter -L$(KLFITTERDIR)/lib -lKLFitterExtras
ifneq ($(MACHINE), $(MAC))
	KLFITTERLIBSO = $(KLFITTERDIR)/lib/libKLFitter.so
	KLFITTERLIBSO_EXTRAS = $(KLFITTERDIR)/lib/libKLFitterExtras.so
else
	KLFITTERLIBSO = $(KLFITTERDIR)/lib/libKLFitter.dylib
	KLFITTERLIBSO_EXTRAS = $(KLFITTERDIR)/lib/libKLFitterExtras.dylib
endif

# Programs
CXX  = g++
RM   = rm -f
ECHO = echo
LN   = ln -s

# flags and libs
CXXFLAGS += $(ROOTCFLAGS) $(BATCFLAGS) -I$(KLFITTERDIR) -I$(KLFITTERDIR_EXTRAS)
CXXFLAGS += -Wall -Wno-deprecated -O2 -ggdb -gstabs
ifneq ($(MACHINE), $(MAC))
	CXXFLAGS += -fPIC
endif
LIBS = $(ROOTLIBS) $(BATLIBS) $(KLFITTERLIBS)

# files
CXXSRCS = runKLFitter_ttH.cxx
CXXOBJS = runKLFitter_ttH.exe
ifeq ($(MACHINE), $(MAC))
	LIBLINK = ../lib/libKLFitter.dylib
	LIBLINK_EXTRAS = ../lib/libKLFitterExtras.dylib
endif
GARBAGE = $(CXXOBJS) $(LIBLINK) $(LIBLINK_EXTRAS)

$(CXXOBJS) : $(CXXSRCS) $(KLFITTERLIBSO) $(KLFITTERLIBSO_EXTRAS)
	@if [ "$(MACHINE)" == "$(MAC)" ]; then if [ ! -e $(LIBLINK) ]; then $(LN) $(KLFITTERLIBSO) $(LIBLINK) ; fi; fi
	@if [ "$(MACHINE)" == "$(MAC)" ]; then if [ ! -e $(LIBLINK_EXTRAS) ]; then $(LN) $(KLFITTERLIBSO_EXTRAS) $(LIBLINK_EXTRAS) ; fi; fi
	$(CXX) $(CXXFLAGS) $(LIBS) $< -o $@

$(KLFITTERLIBSO) : $(KLFITTERLIBSO_EXTRAS)

$(KLFITTERLIBSO_EXTRAS): update-lib

update-lib:
	@cd $(KLFITTERDIR)/scripts && $(MAKE) extras

clean :
	$(RM) $(GARBAGE)
