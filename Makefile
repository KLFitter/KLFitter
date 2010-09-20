# check if OS is linux (Linux) or mac (Darwin)
MACHINE = $(shell uname -s)
LINUX   = Linux
MAC     = Darwin

# ROOT
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs) -lMinuit

# BAT
BATCFLAGS = -I$(BATINSTALLDIR)
BATLIBS   = -L$(BATINSTALLDIR)/lib -lBAT

# programs
CXX    = g++
RM     = rm -f
MKDIR  = mkdir -p
ECHO   = echo

# directories
INCLUDEDIR = include
SRCDIR     = src
OBJDIR     = obj

# flags and libs
CXXFLAGS += $(ROOTCFLAGS) $(BATCFLAGS) -I$(INCLUDEDIR)
CXXFLAGS += -Wall -Wno-deprecated -O2 -ggdb
ifneq ($(MACHINE), $(MAC))
	CXXFLAGS += -fPIC
endif
LIBS = $(ROOTLIBS) $(BATLIBS)
ifneq ($(MACHINE), $(MAC))
	SOFLAGS = -shared     # for linux
else
	SOFLAGS = -dynamiclib # for mac
endif


# files
CXSRCS = ResolutionBase.cxx \
	ResGauss.cxx \
	ResDoubleGaussE.cxx \
	ResDoubleGaussPt.cxx \
	Fitter.cxx \
	DetectorBase.cxx \
	Detector.cxx \
	DetectorAtlas.cxx \
	DetectorDummy.cxx \
	PhysicsConstants.cxx \
	Particles.cxx \
	Permutations.cxx \
	InterfaceBase.cxx \
	InterfaceRoot.cxx \
	InterfaceOutput.cxx \
	InterfaceGoTopTree.cxx \
	InterfaceDummy.cxx \
	LikelihoodBase.cxx \
	LikelihoodTTGamma.cxx \
	LikelihoodTTGamma_RadTopProd.cxx \
	LikelihoodTTGamma_HadTopRadDecay.cxx \
	LikelihoodTTGamma_LepTopRadDecay.cxx \
	LikelihoodTTGamma_HadWRadDecay.cxx \
	LikelihoodTTGamma_LepWRadDecay.cxx \
	LikelihoodTopLeptonJets.cxx \
	LikelihoodTTHElectron.cxx \
	MatchingTool.cxx \
	SelectionTool.cxx \
	PhotonType.cxx \
	ReadConfigFile.cxx

PREPROCESSOR_H = $(INCLUDEDIR)/PREPROC.h
HEADERSRCS     = $(patsubst %.cxx, $(INCLUDEDIR)/%.h, $(CXSRCS))
CXXSRCS        = $(patsubst %.cxx, $(SRCDIR)/%.cxx,   $(CXSRCS))
CXXOBJS        = $(patsubst %.cxx, $(OBJDIR)/%.o,     $(CXSRCS))
ifneq ($(MACHINE), $(MAC))
	LIBSO = libKLFitter.so
else
	LIBSO = libKLFitter.dylib
endif
GARBAGE = $(CXXOBJS) $(LIBSO) $(DYLIB)

# rule for shared library
$(LIBSO) : $(CXXOBJS)
	@if [ "$(MACHINE)" != "$(MAC)" ]; then $(CXX) $(SOFLAGS) $^ -o $(LIBSO); else $(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) $(CXXOBJS) -o $(LIBSO); fi

# rule for single class output files
$(OBJDIR)/%.o : $(SRCDIR)/%.cxx $(INCLUDEDIR)/%.h $(PREPROCESSOR_H)
	@if [ ! -e $(OBJDIR) ]; then $(MKDIR) $(OBJDIR); fi
	$(CXX) $(CXXFLAGS) -c $< -o $@

# rule for making tar ball
tarball : 
	make clean
	@tar -cf KLFitter.tar \
	$(CXXSRCS) \
	$(HEADERSRCS) \
  include/PREPROC.h \
	doc/README_v1.3 \
	examples/top_ejets/setup.sh \
	examples/top_ejets/Makefile \
	examples/top_ejets/runKLFitter.c \
	examples/top_ejets/input.root \
	Makefile 
	@gzip KLFitter.tar 

# rule for making clean
clean :
	$(RM) $(GARBAGE)

# rule for programs and option printing
print :
	$(ECHO) compiler             : $(CXX)
	$(ECHO) c++ srcs             : $(CXXSRCS)
	$(ECHO) pre-processor header : $(PREPROCESSOR_H)
	$(ECHO) c++ objs             : $(CXXOBJS)
	$(ECHO) shared library       : $(LIBSO)
	$(ECHO) c++ flags            : $(CXXFLAGS)
	$(ECHO) libs                 : $(LIBS)
	$(ECHO) so flags             : $(SOFLAGS)
	$(ECHO) machine              : $(MACHINE)
	$(ECHO) linux                : $(LINUX)
	$(ECHO) mac                  : $(MAC)

