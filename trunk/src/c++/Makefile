VERSION    = 1
PATCHLEVEL = 0
SUBLEVEL   = 0
VERBOSE    = 3
BUILDOPTS  =

SAGEROOT=..
include $(SAGEROOT)/config/Global.make

# Common targets are independent of mped/mlocus issues and are common in most SAGE
# applications, listed in dependency order

# Each entry listed in COMPILE order, not LINK order:
#

CORE          = globals containers error LSF fortran numerics
TOOLS         = util output app functors
PEDIGREE_DATA = mped mlocus rped pairs
DATA_CLEANING = fped func sampling gelim data
PEDIGREE_ALGS = peeling lvec mcmc ibd palbase
ALGORITHMS    = maxfun
APIS          = maxfunapi mfsubmodels
APPLICATIONS  = ageon assoc decipher fcor freq genibd lodlink lodpal \
                markerinfo mlod pedinfo relpal reltest segreg sibpal tdtex

TARGETS = $(CORE) $(TOOLS) $(PEDIGREE_DATA) $(DATA_CLEANING) $(PEDIGREE_ALGS) $(ALGORITHMS) $(APIS) $(APPLICATIONS)

#
# Global rules for INTERNAL build nodes.

#
all: subs

subs:
	@unset MAKEFLAGS; set -e; DIR=`pwd` ; for i in $(TARGETS); do \
	  cd $$i && \
	  $(MAKE) $(BUILDOPTS) BUILD=$(BUILD) VERBOSE=$(VERBOSE) && \
	  cd $$DIR; \
	done

install: dummy
	@$(MAKE) subs BUILDOPTS=install BUILD=$(BUILD) VERBOSE=1

depend: dummy
	@$(MAKE) subs BUILDOPTS=depend BUILD=$(BUILD)

test: dummy
	@$(MAKE) subs BUILDOPTS=test BUILD=$(BUILD) VERBOSE=1

testbuild: dummy
	@$(MAKE) subs BUILDOPTS=testbuild BUILD=$(BUILD) VERBOSE=1

docs: dummy
	@echo "[lib categories]"                 >  docs/~~docdict
	@echo "CORE:          $(CORE)"           >> docs/~~docdict
	@echo "TOOLS:         $(TOOLS)"          >> docs/~~docdict
	@echo "PEDIGREE_DATA: $(PEDIGREE_DATA)"  >> docs/~~docdict
	@echo "DATA_CLEANING: $(DATA_CLEANING)"  >> docs/~~docdict
	@echo "PEDIGREE_ALGS: $(PEDIGREE_ALGS)"  >> docs/~~docdict
	@echo "ALGORITHMS:    $(ALGORITHMS)"     >> docs/~~docdict
	@echo "APIS:          $(APIS)"           >> docs/~~docdict
	@echo "APPLICATIONS:  $(APPLICATIONS)"   >> docs/~~docdict

	@../../util/doc/build_main_pages.py

	@rm docs/~~docdict

	@../../util/doc/moddocs.py

	@../../util/doc/create_pdf_manual.py

	@../../util/doc/create_highest_level.py

	@$(MAKE) subs BUILDOPTS=docs BUILD=$(BUILD) VERBOSE=1

relink: dummy
	@$(MAKE) subs BUILDOPTS=relink BUILD=$(BUILD)

clean: cleanlibs clean_targets

clean_targets:
	@-$(MAKE) subs BUILDOPTS=clean TARGETS="$(TARGETS)"

cleanlibs: dummy
	@-rm -f lib/lib*.a
	@-rm -f tmp/*.pch

cleanall: cleanlibs
	@-$(MAKE) subs BUILDOPTS=cleanall TARGETS="$(TARGETS)"

clean_source: cleanlibs
	@-$(MAKE) subs BUILDOPTS=clean_source TARGETS="$(TARGETS)"


dummy:

help:
	@cat $(SAGEROOT)/config/Makefile.help
