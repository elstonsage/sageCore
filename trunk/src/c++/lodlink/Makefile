# Makefile for Multipedigree Stats Objects
#
#  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Apr  5 1996
#             1.01  kbj  Refinements and confusion     Jun 10 1996
#             1.1   kbj  Finally got it right          Jul  1 1996
#
#  Copyright (c) 1996  R.C. Elston
#

#==========================================================================
# General target options                                                  |
#--------------------------------------------------------------------------
#
# -- These options apply to the entire Makefile, not any single target
#    o TARGET  is used for making a specific target
#    o TARGETS is the list of targets to be make by the recursive rule
#    o VERSION and TARPREFIX set parts of the archive name when the
#      tar rule is called.

  LSFROOT=..
  SAGEROOT=../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = "Lodlink"
  TARGET      =
  TARGETS     = liblodlink.a lodlink$(EXE)
  TESTTARGETS = liblodlink.a lodlink$(EXE) test_parser$(EXE) test_mle$(EXE) test_tcalc$(EXE) test_peeler$(EXE) test_mpcalc$(EXE) test_max$(EXE)
  TARPREFIX   = LODLINK
  TESTS       = runall lodlink

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = mle_sub_model.h instructions.h parser.h input.h \
                trans_calculator.h peeler.h likelihood.h max_opt.h \
                linkage_results.h linkage_tests.h analysis.h \
                homogeneity_results.h homogeneity_tests.h definitions.h \
                lod_results.h lods.h output.h geno_elim.h \
                results.h tasks.h genotypes.h genotype_results.h

  SRCS        = mle_sub_model.cpp instructions.cpp parser.cpp input.cpp \
                trans_calculator.cpp peeler.cpp likelihood.cpp max_opt.cpp \
                linkage_results.cpp linkage_tests.cpp analysis.cpp \
                homogeneity_results.cpp homogeneity_tests.cpp definitions.cpp\
                lods.cpp output.cpp geno_elim.cpp \
                results.cpp tasks.cpp genotypes.cpp genotype_results.cpp

  DEP_SRCS    = lodlink.cpp test_parser.cpp test_mle.cpp test_tcalc.cpp \
                test_peeler.cpp test_mpcalc.cpp test_max.cpp

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS = ${SRCS:.cpp=.h,}

    #======================================================================
    #   Target: liblodlink.a                                              |
    #----------------------------------------------------------------------

       liblodlink.a.NAME     = "lodlink library"
       liblodlink.a.TYPE     = LIB
       liblodlink.a.CXXFLAGS = -I../mped
       liblodlink.a.OBJS     = $(OBJS)
       liblodlink.a.CP       = ../lib/liblodlink.a


    #======================================================================
    #   Target: lodlink
    #----------------------------------------------------------------------

       lodlink$(EXE).NAME     = LODLINK Application
       lodlink$(EXE).INSTALL  = yes
       lodlink$(EXE).TYPE     = C++
       lodlink$(EXE).OBJS     = lodlink.o
       lodlink$(EXE).DEP      = liblodlink.a
       lodlink$(EXE).LDLIBS   = -llodlink $(LIB_ALL) 

    #======================================================================
    #   Target: test_parser                                                |
    #----------------------------------------------------------------------

       test_parser$(EXE).NAME     = Test of lodlink parser
       test_parser$(EXE).TYPE     = C++
       test_parser$(EXE).OBJS     = test_parser.o
       test_parser$(EXE).DEP      = liblodlink.a
       test_parser$(EXE).LDLIBS   = -llodlink  $(LIB_ALL) 

    #======================================================================
    #   Target: test_mle                                                  |
    #----------------------------------------------------------------------

       test_mle$(EXE).NAME     = Test of lodlink mle_sub_model
       test_mle$(EXE).TYPE     = C++
       test_mle$(EXE).OBJS     = test_mle.o
       test_mle$(EXE).DEP      = liblodlink.a
       test_mle$(EXE).LDLIBS   = -llodlink $(LIB_ALL) 

    #======================================================================
    #   Target: test_tcalc                                                  |
    #----------------------------------------------------------------------

       test_tcalc$(EXE).NAME     = Test of lodlink trans_calculator class
       test_tcalc$(EXE).TYPE     = C++
       test_tcalc$(EXE).OBJS     = test_tcalc.o
       test_tcalc$(EXE).DEP      = liblodlink.a
       test_tcalc$(EXE).LDLIBS   = -llodlink $(LIB_ALL) 

    #======================================================================
    #   Target: test_peeler                                                  |
    #----------------------------------------------------------------------

       test_peeler$(EXE).NAME     = Test of lodlink peeler class
       test_peeler$(EXE).TYPE     = C++
       test_peeler$(EXE).OBJS     = test_peeler.o
       test_peeler$(EXE).DEP      = liblodlink.a
       test_peeler$(EXE).LDLIBS   = -llodlink $(LIB_ALL) 

    #======================================================================
    #   Target: test_mpcalc                                                  |
    #----------------------------------------------------------------------

       test_mpcalc$(EXE).NAME     = Test of likelihood calculator classes
       test_mpcalc$(EXE).TYPE     = C++
       test_mpcalc$(EXE).OBJS     = test_mpcalc.o
       test_mpcalc$(EXE).DEP      = liblodlink.a
       test_mpcalc$(EXE).LDLIBS   = -llodlink  $(LIB_ALL) 

    #======================================================================
    #   Target: test_max                                                  |
    #----------------------------------------------------------------------

       test_max$(EXE).NAME     = Test of maximization
       test_max$(EXE).TYPE     = C++
       test_max$(EXE).OBJS     = test_max.o
       test_max$(EXE).DEP      = liblodlink.a
       test_max$(EXE).LDLIBS   = -llodlink  $(LIB_ALL) 

    #======================================================================
    #   Target: scratch                                                  |
    #----------------------------------------------------------------------

       scratch$(EXE).NAME     = scratch
       scratch$(EXE).TYPE     = C++
       scratch$(EXE).OBJS     = scratch.o
       scratch$(EXE).DEP      = liblodlink.a
       scratch$(EXE).LDLIBS   = -llodlink  $(LIB_ALL) 

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

