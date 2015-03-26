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

  TARGET_NAME = "Assoc"
  TARGET      =
  TARGETS     = libassoc.a assoc$(EXE)
  TESTTARGETS = libassoc.a assoc$(EXE)
  TARPREFIX   = ASSOC
  TESTS       = runall assoc

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  DATA         = Configuration.cpp Residuals.cpp
  KERNEL       = MatrixDefs.cpp MemberCovariateCalculator.cpp Calcs.cpp
  ANALYSIS_CTL = MaximizationWrapper.cpp AnalysisResults.cpp Calculator.cpp OutputFormatter.cpp
  APPLICATION  = Parser.cpp AppData.cpp

  SRCS         = ${DATA} ${KERNEL} ${ANALYSIS_CTL} ${APPLICATION}

  DEP_SRCS     = assoc.cpp

  OBJS         = ${SRCS:.cpp=.o}

  HEADERS = ${SRCS:.cpp=.h,}

    #======================================================================
    #   Target: libssoc.a                                                |
    #----------------------------------------------------------------------

       libassoc.a.NAME     = "Assoc library"
       libassoc.a.TYPE     = LIB
       libassoc.a.OBJS     = $(OBJS)
       libassoc.a.CXXFLAGS = # -DAUTOTRACE_ON
       libassoc.a.CP       = ../lib/libassoc.a

    #======================================================================  
    #   Target: assoc
    #----------------------------------------------------------------------
                        
       assoc$(EXE).NAME          = Assoc Application
       assoc$(EXE).INSTALL       = yes
       assoc$(EXE).TYPE          = C++
       assoc$(EXE).OBJS          = assoc.o
       assoc$(EXE).DEP           = libassoc.a
       assoc$(EXE).LDLIBS        = -lassoc $(LIB_ALL) 

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

