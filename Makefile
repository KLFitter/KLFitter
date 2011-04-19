MACHINE = $(shell uname -s)
LINUX   = Linux
MAC     = Darwin

LN = ln -s
RM = rm -rf

ifneq ($(MACHINE), $(MAC))
	LIB = libKLFitter.so
	LIB_EXTRAS = libKLFitterExtras.so
else
	LIB = libKLFitter.dylib
	LIB_EXTRAS = libKLFitterExtras.dylib
endif
LIBDIR = library
LIBDIR_EXTRAS = extras

GARBAGE = $(LIB) $(LIB_EXTRAS)

library: $(LIB)

extras: $(LIB) $(LIB_EXTRAS)

$(LIB): $(LIBDIR)/$(LIB)
	$(RM) $@
	$(LN) $<

$(LIB_EXTRAS): $(LIBDIR_EXTRAS)/$(LIB_EXTRAS) $(LIB) 
	$(RM) $@
	$(LN) $<

$(LIBDIR)/$(LIB): update-lib

$(LIBDIR_EXTRAS)/$(LIB_EXTRAS): update-extras

update-lib:
	@cd $(LIBDIR) && $(MAKE)

update-extras:
	@cd $(LIBDIR_EXTRAS) && $(MAKE)

clean:
	$(RM) $(GARBAGE)
	cd $(LIBDIR) && $(MAKE) clean
	cd $(LIBDIR_EXTRAS) && $(MAKE) clean
