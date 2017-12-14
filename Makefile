MACHINE = $(shell uname -s)
LINUX   = Linux
MAC     = Darwin

INCDIR = KLFitter
SRCDIR = Root
OBJDIR = lib
DESTDIR = dest-tmp

CXX = g++
MKDIR = mkdir -p
RM = rm -f
CP = cp -r

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs) -lMinuit

BATCFLAGS = -I$(BATINSTALLDIR)/include
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

CXXFLAGS = $(ROOTCFLAGS) $(BATCFLAGS) -I$/$(INCDIR) -Wall -Wno-deprecated -O2 -ggdb -g
ifneq ($(MACHINE), $(MAC))
	CXXFLAGS += -fPIC
endif
LIBS     = $(ROOTLIBS) $(BATLIBS)

# rule for shared library
$(LIBSO): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) $+ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx $(INCDIR)/%.h
	@if [ ! -e $(OBJDIR) ]; then $(MKDIR) $(OBJDIR); fi
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: all

all: $(LIBSO)

.PHONY: clean

clean:
	$(RM) $(OBJ) $(LIBSO)

.PHONY: install

install: all
	$(MKDIR) $(DESTDIR)/include
	$(MKDIR) $(DESTDIR)/lib
	$(CP) $(LIBSO) $(DESTDIR)/lib/
	$(CP) $(INCDIR) $(DESTDIR)/include/
