# Makefile for New Fcor
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

  TARGET_NAME = MultiPedigree Correlations with Standard Errors
  TARGET      =
  TARGETS     = libfcor.a fcor$(EXE)
  TESTTARGETS = libfcor.a fcor$(EXE)
  TARPREFIX   = FCOR
  TESTS       = runall fcor   #true

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = parser.h           input.h           pairsetcountinfo.h   \
                setweight.h        pairsetdata.h     optimal_weight.h     \
                corcal.h           subcalM.h         accal.h              \
                htest.h            output.h          analysis.h           \
                fcor.h

  SRCS        = parser.cpp         input.cpp         pairsetcountinfo.cpp \
                setweight.cpp      pairsetdata.cpp   optimal_weight.cpp   \
                corcal.cpp         subcalM.cpp       accal.cpp            \
                htest.cpp          output.cpp        analysis.cpp

  DEP_SRCS    = fcor.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libfcor.a                                                 |
    #----------------------------------------------------------------------

       libfcor.a.NAME     = New Fcor Library
       libfcor.a.TYPE     = LIB
       libfcor.a.OBJS     = $(OBJS)
       libfcor.a.CP       = ../lib/libfcor.a

    #======================================================================
    #   Target: fcor                                                      |
    #----------------------------------------------------------------------

       fcor$(EXE).NAME     = "Multi-pedigree correlation program"
       fcor$(EXE).TYPE     = C++
       fcor$(EXE).INSTALL  = yes
       fcor$(EXE).CXXFLAGS = 
       fcor$(EXE).OBJS     = fcor.o
       fcor$(EXE).LDLIBS   = -lfcor $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

