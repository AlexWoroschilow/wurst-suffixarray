# 20 march 2001
# rcsid = $Id: Makefile,v 1.45 2007/06/01 12:22:18 mmosisch Exp $
.POSIX:

LIBNAME=wurst

ECHO=echo

CFLAGS=$(OPTIMIZE) $(CWARN) $(DEFINES) $(CFLAGS_UPPER)

SHELL = /bin/sh

SRCS1 = align.c altscores.c amino_a.c binver.c cmp_dmat.c common.c
SRCS2 = classifyStructure.c class_model.c lsqf.c multialign.c
SRCS3 = coord.c coordinfo.c dbg.c dihedral.c dipep.c e_malloc.c fio.c srand.c
SRCS4 = geo_gap.c matrix.c mgc_num.c misc.c model.c mprintf.c pair_set.c
SRCS5 = pair_set_p.c param_fx.c pair_set_chim.c pdbin.c pdbout.c prob_vec.c
SRCS6 = read_ac_strct.c read_ac.c read_blst.c read_mat.c read_sec.c read_seq.c rescore.c
SRCS7 = scor_set.c score_fx.c score_mat.c score_probvec.c score_sec.c
SRCS8 = score_smat.c scratch.c sec_s.c silly.c str.c seq.c vec.c
#SRCS9 = srand.c

SRCS  = $(SRCS1) $(SRCS2) $(SRCS3) $(SRCS4) $(SRCS5) $(SRCS6) $(SRCS7) $(SRCS8)

WURSTLIB=libwurst.a

nothing:
	echo 'Do not run make from here. cd .. ; make'

all: Makefile.$(LIBNAME)
	$(MAKE) -f Makefile.$(LIBNAME) $(WURSTLIB)
	$(MAKE) -f Makefile.$(LIBNAME) perlstuff

# Next few lines are necessary on shitbox compaq - not on linux
# or solaris.  Can't tell if this is a bug in their make or a
# non-standard feature of sun, linux, whatever.
#COMPILE.c=$(CC) $(CFLAGS) $(CPPFLAGS) -c
#.c.a:
#	$(COMPILE.c) -o $% $<
#	$(AR) $(ARFLAGS) $@ $%
#	$(RM) $%

$(WURSTLIB): gsldir/libminigsl.a $(LIB_DEPS)

Makefile.$(LIBNAME): depend Makefile
	@$(ECHO) "Rebuilding Makefile.$(LIBNAME)"
	@$(ECHO) ".POSIX:"                           >  Makefile.$(LIBNAME)
	@$(ECHO) "# Temporary deletable makefile"    >> Makefile.$(LIBNAME)
	@$(ECHO) "# Do not bother editing this file" >> Makefile.$(LIBNAME)
	@$(ECHO) "# Your changes will go away"       >> Makefile.$(LIBNAME)
	@$(ECHO) "OBJS= " | tr -d \\\012\            >> Makefile.$(LIBNAME)
	@for i in $(SRCS); do \
            $(ECHO)  $$i  | sed -e s/\\\(.*\\\)\.c/\\1.o\ / | tr -d \\\012 \
            >> Makefile.$(LIBNAME);\
        done;\
        $(ECHO)                                     >> Makefile.$(LIBNAME) ;\
        $(ECHO) "LIB_DEPS = " | tr -d \\\012\       >> Makefile.$(LIBNAME)
	@for i in $(SRCS); do \
            $(ECHO) $$i | sed -e s/\\\(.*\\.\\\)c/$(WURSTLIB)\(\\1o\)\ / \
            | tr -d \\\012\
                                                    >> Makefile.$(LIBNAME);\
        done;\
	$(ECHO)                                     >> Makefile.$(LIBNAME);\
	cat depend Makefile                         >> Makefile.$(LIBNAME)

install:
	make all
	(cd Wurst; make $@)

depend: Makefile
	$(DEPENDER) $(SRCS) > depend

perlstuff: Wurst/Makefile
	@cd Wurst; CC=$(CC) OPTIMIZE='$(OPTIMIZE)' $(MAKE)

Wurst/Makefile: Wurst/Makefile.PL
	cd Wurst; CC=$(CC) OPTIMIZE='$(OPTIMIZE)' MANDEST=$(MANDEST)\
                  perl Makefile.PL PREFIX=$(PPREFIX) LIB=$(PLIB)
lint:
	lint -errhdr=%all -Nlevel=4 $(SRCS) -lm

gsldir/libminigsl.a:
	cd gsldir; $(MAKE) DEPENDER='$(DEPENDER)'

splint:
	splint -preproc +posixlib -realcompare $(SRCS)

tags:
	etags $(SRCS)

clean:
	@-rm -f *.o *~ *.ln
	-rm -f Makefile.$(LIBNAME)
	-rm -f $(WURSTLIB)
	@-rm -f *.bbg *.da   # gnu coverage files
	-rm -f depend core TAGS
	@(cd Wurst; $(MAKE) -f Makefile.remedial clean);
	@(cd Wurst; if [ -f Makefile ] ; then $(MAKE) clean; fi )
	@(cd Wurst; if [ -f Makefile ] ; then $(MAKE) realclean; fi )
	@(cd gsldir ; $(MAKE) clean);
