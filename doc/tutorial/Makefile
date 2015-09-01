# $Id: Makefile 15006 2007-08-21 17:29:21Z bangerth $
# Create proper HTML output from the different example programs

# include paths and global variables

#include $(shell echo ${DEALHOME})/common/Make.global_options
PERL = /usr/bin/perl
examples_dir = ../../testsite/

# plain names of those files that we shall treat




example-names                                                                                                                              = step-27 

 
# step-4 step-5 step-6 step-7 step-8 step-9 step-10 step-11 step-12 step-13 step-14 step-15 step-17 \
#  step-21 step-22 step-23 step-24 step-25 step-26 step-27 step-28
#$(notdir $(shell echo $(examples_dir)/step-?)) $(notdir $(shell echo $(examples_dir)/step-??))

# prefix names for generated intermediate files
gen-example-names = $(addprefix generated/, $(example-names))

# names of commented program files, plain files, and tocs for the commented
# and commented and plain versions of the program as well as a combined one
example-toc-prog   = $(addsuffix _toc.prog,$(gen-example-names))
example-toc-plain  = $(addsuffix _toc.plain,$(gen-example-names))
example-toc-intro  = $(addsuffix _toc.intro,$(gen-example-names))
example-toc-results= $(addsuffix _toc.results,$(gen-example-names))
example-toc        = $(addsuffix _toc.combined,$(gen-example-names))

# same for intermediate files for doxygen output
example-dox-prog   = $(addsuffix _prog.dox,$(gen-example-names))
example-dox-plain  = $(addsuffix _plain.dox,$(gen-example-names))


# finally a target for combined Doxygen files
example-doxygen     = $(addprefix doxygen/,$(addsuffix .h,$(example-names)))


# a makefile command that extracts the base name of a program
get_basename = $(shell echo $(1) | $(PERL) -pi -e 's/(generated|doxygen).//g; s/\..*//g; s/_(prog|plain|toc.*)//g;')


all: directories $(example-doxygen)


$(example-toc-prog):
	@echo ================== Making $@
	@cat $(examples_dir)/$(call get_basename, $@)/*.cc \
	  | $(PERL) program2toc \
	  > $@


$(example-toc-plain):
	@echo ================== Making $@
	@cat generated/$(call get_basename, $@)_toc.prog \
	  | $(PERL) -pi -e 's/href=\"\#/href=\"\#plain-/g;' \
	  > $@


$(example-toc-intro):
	@echo ================== Making $@
	@cat $(examples_dir)/$(call get_basename, $@)/doc/intro.dox \
	  | egrep -i "<h[3456]>.*</h[3456]>" \
	  | $(PERL) ./intro2toc \
	  > $@


$(example-toc-results):
	@echo ================== Making $@
	@cat $(examples_dir)/$(call get_basename, $@)/doc/results.dox \
	  | egrep -i "<h[3456]>.*</h[3456]>" \
	  | $(PERL) ./intro2toc \
	  > $@


$(example-toc):
	@echo ================== Making $@
	@echo "<table class=\"tutorial\" width=\"100%\">" > $@
	@echo "<tr><th colspan=\"2\"><b><small>Table of contents</small></b></th></tr>" >> $@
	@echo "<tr><td width=\"50%\" valign=\"top\">" >> $@
	@echo "<ol>" >> $@
	@echo "  <li> <a href=\"#Intro\" class=bold>Introduction</a>" >> $@
	@cat generated/$(call get_basename, $@)_toc.intro >> $@
	@echo "  <li> <a href=\"#CommProg\" class=bold>The commented program</a>" >> $@
	@cat generated/$(call get_basename, $@)_toc.prog >> $@
	@echo "</ol></td><td width=\"50%\" valign=\"top\"><ol>" >> $@
	@echo "  <li value=\"3\"> <a href=\"#Results\" class=bold>Results</a>" >> $@
	@cat generated/$(call get_basename, $@)_toc.results >> $@
	@echo "  <li> <a href=\"#PlainProg\" class=bold>The plain program</a>" >> $@
	@cat generated/$(call get_basename, $@)_toc.plain >> $@
	@echo "</ol> </td> </tr> </table>" >> $@



$(example-dox-prog):
	@echo ================== Making $@
	python gather_cu_and_cpp.py $(examples_dir)/$(call get_basename, $@)/
	@cat $(examples_dir)/$(call get_basename, $@)/step-cu.cc \
	  | $(PERL) program2doxygen \
	  > $@


$(example-dox-plain):
	@echo ================== Making $@
	@echo " * <a name=\"PlainProg\"></a>" > $@
	@echo " * <h1> The plain program</h1>" >> $@
	@echo " * " >> $@
	@echo " * (If you are looking at a locally installed CUDA HPC Praktikum version, then the" >> $@
	@echo " * program can be found at <i>" >> $@
	@echo " * " $(examples_dir)/$(call get_basename,$@)/*.cc | $(PERL) -pi -e 's!/! /!g;' >> $@
	@echo " * </i>. Otherwise, this is only" >> $@
	@echo " * the path on some remote server.)" >> $@
	@echo " @code" >> $@
	@cat $(examples_dir)/$(call get_basename, $@)/*.cc | $(PERL) program2doxyplain >> $@
	@echo " @endcode" >> $@


$(example-doxygen):
	@echo ================== Assembling $@
	@echo "/**" > $@
	@echo " * @page $(shell echo $(call get_basename,$@) | $(PERL) -pi -e 's/-/_/g;') \
	                The $(call get_basename,$@) tutorial program" >> $@
	@echo "@htmlonly" >> $@
	@cat generated/$(call get_basename, $@)_toc.combined >> $@
	@echo "@endhtmlonly" >> $@
	@cat $(examples_dir)/$(call get_basename, $@)/doc/intro.dox \
		| $(PERL) create_anchors >> $@
	@cat generated/$(call get_basename, $@)_prog.dox >> $@
	@cat $(examples_dir)/$(call get_basename, $@)/doc/results.dox \
		| $(PERL) create_anchors >> $@
	@cat generated/$(call get_basename, $@)_plain.dox >> $@
	@echo " */" >> $@


validate-xrefs:
	@$(PERL) ../../scripts/doc/validate-xrefs.pl \
                $(filter-out head.html foot.html, \
                             $(shell echo *.html))


clean:
	-rm -f $(example-toc) \
		$(example-toc-prog) $(example-toc-plain) \
		$(example-toc-intro) $(example-toc-results) \
		$(example-dox-prog) $(example-dox-plain) $(example-doxygen)
	-rm -f Makefile.dep



directories:
	( mkdir -p generated )
	( mkdir -p doxygen ) 


.PHONY: validate-xrefs clean


Makefile.dep: $(shell echo $(examples_dir)/*/) Makefile
	@echo ================== Generating $@
	@for i in $(example-names) ; do \
		echo generated/$${i}_prog.html generated/$${i}_plain.html \
		     generated/$${i}_toc.prog \
		     generated/$${i}_prog.dox generated/$${i}_plain.dox : \
			`echo $(examples_dir)/$${i}/` ; \
		echo generated/$${i}_toc.plain : generated/$${i}_toc.prog ; \
		echo generated/$${i}_toc.combined : \
			generated/$${i}_toc.intro \
			generated/$${i}_toc.results \
			generated/$${i}_toc.prog \
			generated/$${i}_toc.plain ; \
		echo doxygen/$$i.h : generated/$${i}_prog.dox generated/$${i}_plain.dox \
				generated/$${i}_toc.combined \
				$(examples_dir)/$$i/doc/intro.dox \
				$(examples_dir)/$$i/doc/results.dox \
			 ; \
		echo generated/$${i}_toc.intro : \
				$(examples_dir)/$$i/doc/intro.dox \
			 ; \
		echo generated/$${i}_toc.results : \
				$(examples_dir)/$$i/doc/results.dox \
			 ; \
	 done \
	> $@

include Makefile.dep
