# src/Makefile  for a Perl/XS/C module
# Makefile for Wurst::Guts
# 14th June 2004

.POSIX:

MODULE = Guts

MODULE_DIR = Wurst/Guts

DEFINES =-DDONT_SEE_RCS -Dcrazy_gdb_bug

INCLUDES =-I../../../src/

CFLAGS   = $(INCLUDES) $(OPTIMIZE) $(CWARN) $(DEFINES)


SRCS = coord_deletion.c model_pscmp.c pair_setxch.c scoranlys.c mcigraph_io.c \
	pair_settrim.c pair_setshft.c seqprof_trim.c
ECHO = echo

nothing:
	echo 'Do not run make from here. cd ../ ; make'

all: Makefile.$(MODULE)
	$(MAKE) -f Makefile.$(MODULE) perlstuff


Makefile.$(MODULE): depend Makefile
	@$(ECHO) "Rebuilding Makefile.$(MODULE)"
	@$(ECHO) "# Temporary deletable makefile"    >  Makefile.$(MODULE)
	@$(ECHO) "# Do not bother editing this file" >> Makefile.$(MODULE)
	@$(ECHO) "# Your changes will go away"       >> Makefile.$(MODULE)
	@$(ECHO) "OBJS= " | tr -d \\\012\            >> Makefile.$(MODULE)
	@for i in $(SRCS); do \
            $(ECHO)  $$i  | sed -e s/\\\(.*\\\)\.c/\\1.o\ / | tr -d \\\012 \
            >> Makefile.$(MODULE);\
        done;\
        $(ECHO)                                     >> Makefile.$(MODULE) ;\
        $(ECHO) "MOD_OBJ_DEPS = " | tr -d \\\012\       >> Makefile.$(MODULE)
	@for i in $(SRCS); do \
            $(ECHO) $$i | sed -e s/\\\(.*\\.\\\)c/\\1o\ / \
            | tr -d \\\012\
                                                    >> Makefile.$(MODULE);\
        done;\
	$(ECHO)                                     >> Makefile.$(MODULE);\
	cat depend Makefile                         >> Makefile.$(MODULE)

install:
	make all
	(cd $(MODULE_DIR); make $@)

depend: Makefile
	$(DEPENDER) $(CFLAGS) $(SRCS) > depend

perlstuff: $(MOD_OBJ_DEPS) Makefile.$(MODULE) $(MODULE_DIR)/Makefile 
	cd $(MODULE_DIR); CC=$(CC) OPTIMIZE=$(OPTIMIZE) $(MAKE)

$(MODULE_DIR)/Makefile: $(MODULE_DIR)/Makefile.PL
	cd $(MODULE_DIR); CC=$(CC) OPTIMIZE='$(OPTIMIZE)' MANDEST=$(MANDEST) \
                MOD_OBJS="$(OBJS)" \
                perl Makefile.PL PREFIX=$(PPREFIX)

lint:
	lint -errhdr=%all -Nlevel=4 $(SRCS) -lm

tags:
	etags $(SRCS)

clean:
	-rm -f *~
	-rm -f *.o
	-rm -f *.ln
	-rm -f Makefile.$(MODULE)
	-rm -f depend
	-rm -f core
	-rm -f TAGS
	(cd $(MODULE_DIR); if [ -f Makefile ] ; then echo Cleaning $(MODULE_DIR); \
	$(MAKE) realclean; fi )

