# Wurst::AlDb 
# makefile - derived from the non-XS modeller model


.POSIX:

MODULE_DIR = Wurst

DEFINES =-DDONT_SEE_RCS -Dcrazy_gdb_bug

INCLUDES =-I../src

CFLAGS   = $(INCLUDES) $(OPTIMIZE) $(CWARN) $(DEFINES)
# This next is coded into gnu make.
# No idea if it is standard or not.

LOADLIBES=-lm

OBJECT=


PERLSCRIPTS=

all: $(OBJECT) perlstuff

perlstuff: $(MODULE_DIR)/Makefile
	cd $(MODULE_DIR); CC=$(CC) OPTIMIZE='$(OPTIMIZE)' $(MAKE)

install: perlstuff
	cd $(MODULE_DIR); CC=$(CC) OPTIMIZE='$(OPTIMIZE)' $(MAKE) $@

perlscripts: 
	if [ ! -d ibin ]; then mkdir ibin; fi; \
	perlpref=`echo $(PPREFIX) | sed -e 's/\\//\\\\\\//g'`; \
	for i in $(PERLSCRIPTS); do \
	  sed -e \
 's/\#__wurstli__/use\ lib\ \"'$$perlpref'\/lib\/perl5\/site_perl\"\;\#/' < $$i >ibin/$$i ;\
	done

$(MODULE_DIR)/Makefile: $(MODULE_DIR)/Makefile.PL
	cd $(MODULE_DIR); CC=$(CC) OPTIMIZE='$(OPTIMIZE)' MANDEST=$(MANDEST)\
                  perl Makefile.PL PREFIX=$(PPREFIX)

clean:
	-rm -f *~
	-rm -f *.o
	-rm -f del_me*
	-rm -f *.ln
	-rm -f core
	-rm -f TAGS
	(cd $(MODULE_DIR); if [ -f Makefile ] ; then $(MAKE) realclean; fi )




