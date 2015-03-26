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

  SAGEROOT=../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = "Age of Onset"
  TARGET      =
  TARGETS     = libageon.a ageon$(EXE)
  TESTTARGETS = libageon.a ageon
  TARPREFIX   = AO
  TESTS       = runall ageon

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS        = AnalysisOutput.cpp            \
                AnalysisWrapper.cpp           \
                Calculator.cpp                \
                AppData.cpp                   \
                Kernel.cpp                    \
                MemberCovariateCalculator.cpp \
                ExtraOutput.cpp               \
                Model.cpp                     \
                Parser.cpp                    \
                Validator.cpp

  DEP_SRCS    = ageon.cpp

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS = ${SRCS:.cpp=.h,}

    #======================================================================
    #   Target: libageon.a                                                |
    #----------------------------------------------------------------------

       libageon.a.NAME     = "Age of Onset library"
       libageon.a.TYPE     = LIB
       libageon.a.CXXFLAGS = -I../mped
       libageon.a.OBJS     = $(OBJS)
       libageon.a.CP       = ../lib/libageon.a

    #======================================================================
    #   Target: ageon
    #----------------------------------------------------------------------

       ageon$(EXE).NAME          = Age of Onset Application
       ageon$(EXE).INSTALL       = yes
       ageon$(EXE).TYPE          = C++
       ageon$(EXE).OBJS          = ageon.o
       ageon$(EXE).DEP           = libageon.a
       ageon$(EXE).LDLIBS        = -lageon $(LIB_ALL) 

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

