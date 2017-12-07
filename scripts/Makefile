MACHINE = $(shell uname -s)
LINUX   = Linux
MAC     = Darwin

DIR = ..

INCDIR = $(DIR)/KLFitter
SRCDIR = $(DIR)/Root
OBJDIR = $(DIR)/obj
LIBDIR = $(DIR)/lib

DIR_EXTRA = $(DIR)/KLFitterExtras

INCDIR_EXTRA = $(DIR_EXTRA)/KLFitterExtras
SRCDIR_EXTRA = $(DIR_EXTRA)/Root
OBJDIR_EXTRA = $(DIR_EXTRA)/obj

CXX = g++
RM = rm -f

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lMinuit

BATCFLAGS = -I$(BATINSTALLDIR)
BATLIBS = -L$(BATINSTALLDIR)/lib -lBAT

LIBS = $(ROOTLIBS) $(BATLIBS)

SRC = $(wildcard $(SRCDIR)/*.cxx)
OBJ = $(SRC:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)
SRC_EXTRA = $(wildcard $(SRCDIR_EXTRA)/*.cxx)
OBJ_EXTRA = $(SRC_EXTRA:$(SRCDIR_EXTRA)/%.cxx=$(OBJDIR_EXTRA)/%.o)

CXXFLAGS = $(ROOTCFLAGS) $(BATCFLAGS) -I$(DIR) -Wall -Wno-deprecated -O2 -ggdb -g

ifneq ($(MACHINE), $(MAC))
	LIBSO = $(LIBDIR)/libKLFitter.so
	LIBSO_EXTRA = $(LIBDIR)/libKLFitterExtras.so
	SOFLAGS = -shared
	CXXFLAGS += -fPIC
else
	LIBSO = $(LIBDIR)/libKLFitter.dylib
	LIBSO_EXTRA = $(LIBDIR)/libKLFitterExtras.dylib
	SOFLAGS = -dynamiclib
endif

GARBAGE = $(OBJ) $(LIBSO) $(OBJ_EXTRA) $(LIBSO_EXTRA)

$(LIBSO) : $(OBJ)
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) $+ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx $(INCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

extras : $(LIBSO_EXTRA)

$(LIBSO_EXTRA) : $(OBJ_EXTRA) $(LIBSO) 
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) $+ -o $@

$(OBJDIR_EXTRA)/%.o: $(SRCDIR_EXTRA)/%.cxx $(INCDIR_EXTRA)/%.h
	$(CXX) $(CXXFLAGS) -c -I$(DIR_EXTRA) $< -o $@

clean:
	$(RM) $(GARBAGE)
