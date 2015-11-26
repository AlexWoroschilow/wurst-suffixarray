# 11 April 2000

.POSIX:

SUBDIRS=optim_seq

all:
	@echo Finished in scripts dir
	@for i in $(SUBDIRS); do\
          (cd $$i; echo making $@ in $$i; $(MAKE) $@) ;\
        done

test:
	perl all.pl

clean:
	-rm -f del_me* delete_me.* tmp.*
	@-rm -f *~
	@-rm core
	-rm splx_out*.out
	@for i in $(SUBDIRS); do\
          (cd $$i; echo making $@ in $$i; $(MAKE) $@) ;\
        done
