# 20 march 2001
# rcsid = $Id: Makefile,v 1.2 2006/05/10 09:09:07 torda Exp $
.POSIX:

LIBNAME=minigsl

ECHO=echo

CFLAGS=$(OPTIMIZE) $(CWARN) $(DEFINES) $(CFLAGS_UPPER)

SHELL = /bin/sh

SRCS1 = blas.c error.c
SRCS2 = 
SRCS3 = permutation_init.c
SRCS4 = lu.c
SRCS5 = init.c view.c
SRCS6 = strerror.c block_init.c
SRCS7 = matrix.c matrix_view.c matrix_init.c matrix_rowcol.c matrix_swap.c

SRCS  = $(SRCS1) $(SRCS2) $(SRCS3) $(SRCS4) $(SRCS5) $(SRCS6) $(SRCS7)


WURSTLIB=libminigsl.a

all: Makefile.$(LIBNAME)
	$(MAKE) -f Makefile.$(LIBNAME) $(WURSTLIB)

# Next few lines might be necessary on some machines like compaq.
#COMPILE.c=$(CC) $(CFLAGS) $(CPPFLAGS) -c
#.c.a:
#	$(COMPILE.c) -o $% $<
#	$(AR) $(ARFLAGS) $@ $%
#	$(RM) $%

$(WURSTLIB): $(LIB_DEPS)

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

depend: Makefile
	$(DEPENDER) $(SRCS) > depend

perlstuff: Wurst/Makefile
	cd ../Wurst; CC=$(CC) OPTIMIZE=$(OPTIMIZE) $(MAKE)

tags:
	etags $(SRCS)

splint:
	splint -preproc +posixlib -realcompare $(SRCS)

clean:
	@-rm -f *.o *~ *.ln
	-rm -f Makefile.$(LIBNAME)
	-rm -f $(WURSTLIB)
	@-rm -f *.bbg *.da   # gnu coverage files
	-rm -f depend core TAGS
