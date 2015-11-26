# Perl/XS/C module makefile
# Follows model used in wurst/ and wurst/src
# rcsid = $Id: Makefile,v 1.2 2004/06/29 08:39:26 procter Exp $


#                     PASS-THROUGH MAKEFILE
# ----------------------- Leave alone -------------------------------
# The stuff below should not be tinkered with, unless you know what
# you are doing

SUBDIRS=src

# have to quote because this is a second-level passthrough
FLAGS_TO_PASS=CC=$(CC) OPTIMIZE=$(OPTIMIZE) CWARN=$(CWARN) \
DEFINES=$(DEFINES) PLIB=$(PLIB) PPREFIX=$(PPREFIX) MANDEST=$(MANDEST) \
DEPENDER=$(DEPENDER)

.POSIX:

all:
	@for i in $(SUBDIRS); do\
	  (cd $$i; echo making $@ in $$i; $(MAKE) $@) ;\
	done

install:
	$(MAKE) all;
	(cd src; $(MAKE) $@)

clean:
	@for i in $(SUBDIRS); do\
	  (cd $$i; echo making $@ in $$i; $(MAKE) $@ ) ;\
	done

home:
	$(MAKE) clean;
	( cd ..; tar cXf excl_list - wurst | gzip -c > w.tgz)
