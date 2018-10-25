# Copyright (c) 2009--2018, the KLFitter developer team
#
# This file is part of KLFitter.
#
# KLFitter is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# KLFitter is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
#
INCDIR = include
SRCDIR = src
OBJDIR = obj
LIBDIR = lib
DESTDIR = build

MAINDIR = util
MAINTARGETDIR = bin

TESTDIR = tests
TESTTARGETDIR = test-bin

CXX = g++
MKDIR = mkdir -p
RM = rm -f
CP = cp -r
AR = ar rvs

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs) -lMinuit

BATCFLAGS = -I$(BATINSTALLDIR)/include
BATLIBS   = -L$(BATINSTALLDIR)/lib -lBAT

SRC = $(wildcard $(SRCDIR)/*.cxx)
OBJ = $(SRC:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)
LIBSO = $(LIBDIR)/libKLFitter.so
LIBA = $(LIBDIR)/libKLFitter.a

MAINSRC = $(wildcard $(MAINDIR)/*.cxx)
MAINEXE = $(MAINSRC:$(MAINDIR)/%.cxx=$(MAINTARGETDIR)/%.exe)

TESTSRC = $(wildcard $(TESTDIR)/*.cxx)
TESTEXE = $(TESTSRC:$(TESTDIR)/%.cxx=$(TESTTARGETDIR)/%.exe)

SOFLAGS = -shared
CXXFLAGS = $(ROOTCFLAGS) $(BATCFLAGS) -I$(INCDIR) -Wall -pedantic -O2 -g -std=c++11 -fPIC
LIBS     = $(ROOTLIBS) $(BATLIBS)

# rule for main executables
$(MAINTARGETDIR)/%.exe: $(MAINDIR)/%.cxx $(LIBA)
	$(MKDIR) $(MAINTARGETDIR)
	$(CXX) $(CXXFLAGS) $+ $(LIBS) -o $@

# rule for test executables
$(TESTTARGETDIR)/%.exe: $(TESTDIR)/%.cxx $(LIBA)
	$(MKDIR) $(TESTTARGETDIR)
	$(CXX) $(CXXFLAGS) $+ $(LIBS) -o $@

# rule for shared library
$(LIBSO): $(OBJ)
	@if [ ! -e $(LIBDIR) ]; then $(MKDIR) $(LIBDIR); fi
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) $+ -o $@

# rule for static library
$(LIBA): $(OBJ)
	$(AR) $@ $+

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	@if [ ! -e $(OBJDIR) ]; then $(MKDIR) $(OBJDIR); fi
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: all

all: $(LIBSO) $(LIBA) $(TESTEXE) $(MAINEXE)

.PHONY: examples

examples: $(MAINEXE)

.PHONY: tests

tests: $(TESTEXE)

.PHONY: clean

clean:
	$(RM) $(OBJ) $(LIBSO) $(LIBA)
	$(RM) -r $(MAINTARGETDIR) $(TESTTARGETDIR) $(DESTDIR)
	$(RM) $(TESTEXE)

.PHONY: install

install: all
	$(MKDIR) $(DESTDIR)/include
	$(MKDIR) $(DESTDIR)/lib
	$(MKDIR) $(DESTDIR)/test-bin
	$(CP) $(LIBDIR) $(DESTDIR)/
	$(CP) $(INCDIR) $(DESTDIR)/
	$(CP) $(TESTEXE) $(DESTDIR)/test-bin
