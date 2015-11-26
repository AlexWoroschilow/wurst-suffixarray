# 31 march 2002

.POSIX:

DSSP_DIR=Dssp
DSSPLIB=libdssp.a
DEFINES =-DDONT_SEE_RCS

INCLUDES =-I../src

CFLAGS   = $(INCLUDES) $(OPTIMIZE) $(CWARN) $(DEFINES) $(CFLAGS_UPPER)

OBJECT=dssp.o

$(DSSPLIB): $(OBJECT)
	ar r $(DSSPLIB) $(OBJECT)

all: $(DSSPLIB) perlstuff

perlstuff: $(DSSP_DIR)/Makefile
	cd $(DSSP_DIR); CC=$(CC) OPTIMIZE='$(OPTIMIZE)' $(MAKE)

$(DSSP_DIR)/Makefile: $(DSSP_DIR)/Makefile.PL
	cd $(DSSP_DIR); CC=$(CC) OPTIMIZE='$(OPTIMIZE)' MANDEST=$(MANDEST)\
                  perl Makefile.PL PREFIX=$(PPREFIX)

clean:
	-rm -f *~
	-rm -f *.o
	-rm -f del_me*
	-rm -f *.ln
	-rm -f $(DSSPLIB)
	-rm -f core
	-rm -f TAGS
	(cd Dssp; if [ -f Makefile ] ; then $(MAKE) realclean; fi )




