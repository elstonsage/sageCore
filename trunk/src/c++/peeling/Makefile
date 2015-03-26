# Makefile for Pedigree Peeling
#
#  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Apr  5 1996
#             1.01  kbj  Refinements and confusion     Jun 10 1996
#             1.1   kbj  Finally got it right          Jul  1 1996
#             2.0   gcw  Recreated everything, moved   Dec 30 1997
#
#  Copyright (c) 2001  R.C. Elston
#

#=========================================================================
# TESTS:                                                                 |
#-------------------------------------------------------------------------
#
#   Not all of the genotype elimiation library is tested here.  Additional
#   tests can be found in the relpal/pairs directory.
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
  LSFROOT=..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = "Peeling Library"
  TARGET      =
  TARGETS     = 
  TESTTARGETS = test_peeler2$(EXE)
  TESTS       = 
  VERSION     = 1.0
  TARPREFIX   = PL
  TESTS       = runall peeling

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  DEP_SRCS    = test_peeler2.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: test_peeler2                                                |
    #----------------------------------------------------------------------

       test_peeler2$(EXE).NAME     = Test of MP File Input
       test_peeler2$(EXE).TYPE     = C++
       test_peeler2$(EXE).OBJS     = test_peeler2.o
       test_peeler2$(EXE).LDLIBS   = $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

