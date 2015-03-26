# Makefile for Relationship test Objects
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

  LSFROOT = ..
  SAGEROOT= ../..

  include  $(SAGEROOT)/config/Global.make

  TARGET_NAME = relationship test objects
  TARGET      =
  TARGETS     = libreltest.a reltest$(EXE)
  TESTTARGETS = libreltest.a reltest$(EXE)
  VERSION     = 1.0
  TARPREFIX   = RELTEST
  TESTS       = runall reltest

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = input.h                parser.h        analysis.h     \
                putative_pair.h \
                l2_procedure.h         output.h         reltest.h

  SRCS        = input.cpp              parser.cpp      analysis.cpp     \
                l2_procedure.cpp       output.cpp

  DEP_SRCS    = reltest.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libreltest.a                                              |
    #----------------------------------------------------------------------

       libreltest.a.NAME     = "Relationship Test program Library"
       libreltest.a.TYPE     = LIB
       libreltest.a.CXXFLAGS = #-DDEBUG_VERBOSE -I. -g
       libreltest.a.LDFLAGS  = # --backend -non_shared  #-DDEBUG_VERBOSE -I. -g
       libreltest.a.OBJS     = ${OBJS}
       libreltest.a.CP       = ../lib/libreltest.a

    #======================================================================
    #   Target: reltest                                                   |
    #----------------------------------------------------------------------

       reltest$(EXE).NAME     = "Relationship Test program"
       reltest$(EXE).INSTALL  = yes
       reltest$(EXE).TYPE     = C++
       reltest$(EXE).CXXFLAGS = #-DDEBUG_VERBOSE -I. -g
       reltest$(EXE).LDFLAGS  = # --backend -non_shared  #-DDEBUG_VERBOSE -I. -g
       reltest$(EXE).OBJS     = reltest.o
       reltest$(EXE).LDLIBS   = -lreltest $(LIB_ALL)

include  $(SAGEROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.





-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

