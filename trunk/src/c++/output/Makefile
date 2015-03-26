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

  TARGET_NAME = "Output"
  TARGET      =
  TARGETS     = liboutput.a
  TESTTARGETS = liboutput.a
#  TESTTARGETS = liboutput.a otest$(EXE)
  TARPREFIX   = 
  TESTS       = true
#  TESTS       = runall output

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS        = Table.cpp RenderingRules.cpp Views.cpp
  OBJS        = ${SRCS:.cpp=.o}

  DEP_SRCS    = otest.cpp

  HEADERS     = ${SRCS:.cpp=.h,}

    #======================================================================
    #   Target: liboutput.a                                               |
    #----------------------------------------------------------------------

       liboutput.a.NAME     = "Output library"
       liboutput.a.TYPE     = LIB
       liboutput.a.CXXFLAGS = -I../mped
       liboutput.a.OBJS     = $(OBJS)
       liboutput.a.CP       = ../lib/liboutput.a

    #======================================================================  
    #   Target: otest
    #----------------------------------------------------------------------
                        
       otest$(EXE).NAME          = Output Application
       otest$(EXE).INSTALL       = yes
       otest$(EXE).TYPE          = C++
       otest$(EXE).CXXFLAGS      =
       otest$(EXE).OBJS          = otest.o
       otest$(EXE).DEP           = liboutput.a
       otest$(EXE).LDLIBS        = -loutput $(LIB_EXT_APIS) # $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

