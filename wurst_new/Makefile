# 20 march 2001
# rcsid = $Id: Makefile,v 1.24 2007/05/07 13:31:11 tmargraf Exp $
.POSIX:


# The final perl library will be installed here
# You can use tilde "~", because the perl installation stuff
# knows what to do with it.
PPREFIX    = "~/pl"
PLIB    = "$(PPREFIX)/lib"
MANDEST =$(HOME)/man

DEPENDER='cc -MM'
# Compaq:
#DEPENDER='cc -M'
# Sun / solaris
#DEPENDER = 'cc -xM1'
SHELL = /bin/sh
#CC=cc
#OPTIMIZE=-O2
OPTIMIZE=-g
CFLAGS_UPPER='-fpic -I..'
#CFLAGS_UPPER=-I..
CWARN  ='-pedantic -Wall -Wshadow -Wwrite-strings -W -Wmissing-prototypes -Wno-char-subscripts'
#CWARN  ='-Xc -v'
DEFINES=-DDONT_SEE_RCS
GSL_LIB_FLAGS=-lgsl -lgslcblas 

# ----------------------- Leave alone -------------------------------
# The stuff below should not be tinkered with, unless you know what
# you are doing
SUBDIRS=src/gsldir src scripts pod seq struct dssp 
#EXTRADIRS=extras/Modeller extras/AlDb extras/Pack extras/Guts
EXTRADIRS=
FLAGS_TO_PASS=CC=$(CC) OPTIMIZE=$(OPTIMIZE) CWARN=$(CWARN) \
CFLAGS_UPPER=$(CFLAGS_UPPER) \
DEFINES=$(DEFINES) PLIB=$(PLIB) PPREFIX=$(PPREFIX) MANDEST=$(MANDEST) \
DEPENDER=$(DEPENDER)

all:
	@for i in $(SUBDIRS) $(EXTRADIRS); do\
	  (cd $$i; echo making $@ in $$i; $(FLAGS_TO_PASS) $(MAKE)  $@) ;\
	done

install:
	make all;
	(cd src; make $(FLAGS_TO_PASS) $@)

instxtras:
	make install;
	@for i in $(EXTRADIRS); do\
	  (cd $$i; echo make-install in $$i; $(FLAGS_TO_PASS) $(MAKE) install) ;\
	done


clean:
	@for i in $(SUBDIRS) $(EXTRADIRS); do\
	  (cd $$i; echo making $@ in $$i; $(MAKE) $(FLAGS_TO_PASS) $@) ;\
	done

home:
	make clean
	( cd ..; tar cXf excl_list - wurst | gzip -c > w.tgz)
