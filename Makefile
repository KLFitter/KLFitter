MACHINE = $(shell uname -s)
LINUX   = Linux
MAC     = Darwin

INCDIR = KLFitter
SRCDIR = Root
OBJDIR = lib

CXX = g++
MKDIR = mkdir -p
RM = rm -f

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs) -lMinuit

BATCFLAGS = -I$(BATINSTALLDIR)
BATLIBS   = -L$(BATINSTALLDIR)/lib -lBAT

SRC = $(wildcard $(SRCDIR)/*.cxx)
OBJ = $(SRC:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)
MAIN = $(wildcard *.c)
ifneq ($(MACHINE), $(MAC))
	LIBSO = libKLFitter.so
	SOFLAGS = -shared
else
	LIBSO = libKLFitter.dylib
	SOFLAGS = -dynamiclib
endif

GARBAGE = $(OBJ) $(LIBSO)

CXXFLAGS = $(ROOTCFLAGS) $(BATCFLAGS) -I$/$(INCDIR) -Wall -Wno-deprecated -O2 -ggdb -g
ifneq ($(MACHINE), $(MAC))
	CXXFLAGS += -fPIC
endif
LIBS     = $(ROOTLIBS) $(BATLIBS)

# rule for shared library
$(LIBSO) : $(OBJ)
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) $+ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx $(INCDIR)/%.h
	@if [ ! -e $(OBJDIR) ]; then $(MKDIR) $(OBJDIR); fi
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(GARBAGE)
