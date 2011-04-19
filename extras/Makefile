MACHINE = $(shell uname -s)
LINUX   = Linux
MAC     = Darwin

BASEINCDIR = ../library/include

INCDIR = include
SRCDIR = src
OBJDIR = obj

CXX = g++
MKDIR = mkdir -p
RM = rm -f

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs) -lMinuit

BATCFLAGS = -I$(BATINSTALLDIR)
BATLIBS   = -L$(BATINSTALLDIR)/lib -lBAT
KLFITTERLIBS = -L.. -lKLFitter

SRC = $(wildcard $(SRCDIR)/*.cxx)
OBJ = $(SRC:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o) $(wildcard ../libary/obj/*.o)
MAIN = $(wildcard *.c)
ifneq ($(MACHINE), $(MAC))
	LIBSO = libKLFitterExtras.so
	SOFLAGS = -shared
else
	LIBSO = libKLFitterExtras.dylib
	SOFLAGS = -dynamiclib
endif

GARBAGE = $(OBJ) $(LIBSO)

CXXFLAGS = $(ROOTCFLAGS) $(BATCFLAGS) -I$/$(INCDIR) -I$/$(BASEINCDIR) -Wall -Wno-deprecated -O2 -ggdb -g
ifneq ($(MACHINE), $(MAC))
	CXXFLAGS += -fPIC
endif
LIBS     = $(ROOTLIBS) $(BATLIBS) $(KLFITTERLIBS)

# rule for shared library
$(LIBSO) : $(OBJ)
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) $+ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx $(INCDIR)/%.h
	@if [ ! -e $(OBJDIR) ]; then $(MKDIR) $(OBJDIR); fi
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(GARBAGE)