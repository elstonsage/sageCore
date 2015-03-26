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

  TARGET_NAME = MultiPedigree statistics
  TARGET      = 
  TARGETS     = libpedinfo.a pedinfo$(EXE)
  TESTTARGETS = libpedinfo.a pedinfo$(EXE)
  TARPREFIX   = PEDINFO
  TESTS       = runall pedinfo

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = pedparse.h   pedinfo_input.h   pedinfo.h  stats.h  stats_view.h
  SRCS        = pedparse.cpp pedinfo_input.cpp  stats.cpp   stats_view.cpp
  DEP_SRCS    = pedinfo.cpp
  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libpedinfo.a                                              |
    #----------------------------------------------------------------------

       libpedinfo.a.NAME     = "Pedinfo library"
       libpedinfo.a.TYPE     = LIB
       libpedinfo.a.OBJS     = $(OBJS)
       libpedinfo.a.CP       = ../lib/libpedinfo.a

    #======================================================================
    #   Target: pedinfo                                                   |
    #----------------------------------------------------------------------

       pedinfo$(EXE).NAME     = "Multi-pedigree statistics test program"
       pedinfo$(EXE).INSTALL  = yes
       pedinfo$(EXE).TYPE     = C++
       pedinfo$(EXE).OBJS     = pedinfo.o
       pedinfo$(EXE).LDLIBS   = -lpedinfo $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

