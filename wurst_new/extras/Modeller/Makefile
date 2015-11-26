# Makefile for a module with no XS component
# 31 march 2002

.POSIX:

MODULE_DIR = Wurst/Modeller

DEFINES =-DDONT_SEE_RCS -Dcrazy_gdb_bug

INCLUDES =-I../src

CFLAGS   = $(INCLUDES) $(OPTIMIZE) $(CWARN) $(DEFINES)
# This next is coded into gnu make.
# No idea if it is standard or not.

LOADLIBES=-lm
PERLSCRIPTS=modellermodel.pl
OBJECT=
WURSTINSTALL= `pwd`"/../../"

all: perlstuff

install: perlstuff
	cd $(MODULE_DIR); CC=$(CC) OPTIMIZE='$(OPTIMIZE)' $(MAKE) $@

perlstuff: $(MODULE_DIR)/Makefile perlscripts
	cd $(MODULE_DIR); CC=$(CC) OPTIMIZE='$(OPTIMIZE)' $(MAKE)

# jiffy to insert the exact library path for Wurst into any
# dependent scripts.

perlscripts: 
	@if [ ! -d ibin ]; then mkdir ibin; fi; \
	if [ -n $(PPREFIX) ]; then \
	perlpref=`echo $(PPREFIX) | sed -e 's/\\//\\\\\\//g'`; \
	winstall=`pwd | sed -e 's/\\//\\\\\\//g'`'\/..\/..\/'; \
	for i in $(PERLSCRIPTS); do \
	  sed -e 's/\#__wurstli__/\$$wurst_installpath\ \=\ \"'$$winstall'\"\;\ use\ lib\ \"'$$perlpref'\/lib\/perl5\/site_perl\"\;\#/' < $$i >ibin/$$i ;\
	  echo \*\*README\*\* Relocatable version of $$i in ibin/$$i; \
	done; \
	else echo Lib Paths not updated in ibin/$(PERLSCRIPTS) : PPREFIX not defined. ;\
	fi;

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




