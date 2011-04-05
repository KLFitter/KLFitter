LN = ln -s
RM = rm -rf

LIB = libKLFitter.so
LIBDIR = library

LIB_EXTRAS = libKLFitterExtras.so
LIBDIR_EXTRAS = extras

GARBAGE = $(LIB) $(LIB_EXTRAS)

library: $(LIB)

extras: $(LIB_EXTRAS) $(LIB)

$(LIB): $(LIBDIR)/$(LIB)
	$(RM) $@
	$(LN) $<

$(LIB_EXTRAS): $(LIBDIR_EXTRAS)/$(LIB_EXTRAS) $(LIB) 
	$(RM) $@
	$(LN) $<

$(LIBDIR)/$(LIB): update-lib

$(LIBDIR_EXTRAS)/$(LIB_EXTRAS): update-extras

update-lib:
	cd $(LIBDIR) && $(MAKE)

update-extras:
	cd $(LIBDIR_EXTRAS) && $(MAKE)

clean:
	$(RM) $(GARBAGE)
	cd $(LIBDIR) && $(MAKE) clean
	cd $(LIBDIR_EXTRAS) && $(MAKE) clean
