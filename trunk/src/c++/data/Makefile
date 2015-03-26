# Makefile for Pair Generator Main Body
#
#  Authors: Kevin Jacobs (jacobs@darwin.cwru.edu)
#           Geoff Wedig  (wedig@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Apr  5 1996
#             1.01  kbj  Refinements and confusion     Jun 10 1996
#             1.1   kbj  Finally got it right          Jul  1 1996
#             1.11  gcw  Relative Pair port            Jul 17 1996
#             1.2   gcw  Added Global.make, Rules.make Aug 19 1996
#             1.21  gcw  Written for RPped             Aug 20 1996
#             1.3   gcw  New Source handling           Jan 06 1997
#
#  Copyright (c) 1996  R.C. Elston
#
#  NOTE:  All S.A.G.E. programmers; this Makefile is a good example of how
#         to derive classes from LSF and build them.

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
  include $(SAGEROOT)/config/Global.make

  TTARGETS    = 

  TARGET_NAME = Basic SAGE application data and tools
  TARGET      =
  TARGETS     = libdata.a 
  TESTTARGETS = libdata.a test_cmdline$(EXE)
  VERSION     = 1.0
  TESTS       = runall data

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS = SAGEdata.cpp ArgumentRuleset.cpp

  DEP_SRCS = test_cmdline.cpp

  OBJS = ${SRCS:%.cpp=%.o}

  #======================================================================
  #   Target: libdata.a                                                 |
  #----------------------------------------------------------------------

       libdata.a.NAME     = "SAGE Application Data"
       libdata.a.OBJS     = $(OBJS)
       libdata.a.TYPE     = LIB
       libdata.a.CP       = ../lib/libdata.a

  #----------------------------------------------------------------------
  #   Target: test_cmdline$(EXE)                                        |
  #----------------------------------------------------------------------

       test_cmdline$(EXE).NAME      = Tests the command line parsing stuff
       test_cmdline$(EXE).TYPE      = C++
       test_cmdline$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE
       test_cmdline$(EXE).OBJS      = test_cmdline.o 
       test_cmdline$(EXE).LDFLAGS   = -L.
       test_cmdline$(EXE).DEP       = libdata.a
       test_cmdline$(EXE).LDLIBS    = $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

