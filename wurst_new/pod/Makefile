# 11 Sep 2001
# rcsid $Id: Makefile,v 1.4 2006/02/15 08:18:41 torda Exp $

.POSIX:

SRCS=wurst.pod
MANS=wurst.man coding.man

# Program names  ------------------------------------------
POD2MAN =   pod2man
POD2HTML =  pod2html
POD2LATEX = pod2latex

all: man

install: all
	@echo I should copy the man files to $(MANDEST)

man: $(MANS)

clean:
	-rm -f *.man *.html *~

# Rules only in this section ------------------------------

.SUFFIXES: .man .pod .html

.pod.man:
	$(POD2MAN) $< > $*.man

.pod.html:
	$(POD2HTML) $< --outfile=$*.html
