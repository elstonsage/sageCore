# Makefile for Genotype Elimination
#
#  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Apr  5 1996
#             1.01  kbj  Refinements and confusion     Jun 10 1996
#             1.1   kbj  Finally got it right          Jul  1 1996
#             2.0   gcw  Recreated everything, moved   Dec 30 1997
#
#  Copyright (c) 1996  R.C. Elston
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

  TARGET_NAME = "Genotype Elimination Model"
  TARGET      =
  TARGETS     = libgelim.a
  TESTTARGETS = libgelim.a test_gelim
  TESTS       = runall gelim
  VERSION     = 1.0
  TARPREFIX   = GE

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS        = inconsistency_handler.cpp      \
                valid_parental_genotypes.cpp   \
                geno_eliminate.cpp             \
                ped_imodel_gen.cpp             \
                pedigree_region.cpp

  DEP_SRCS    = test_gelim.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libgelim.a                                                |
    #----------------------------------------------------------------------

       libgelim.a.NAME     = "Genotype Elimination Library"
       libgelim.a.OBJS     = $(OBJS)
       libgelim.a.TYPE     = LIB
       libgelim.a.CP       = ../lib/libgelim.a

    #======================================================================
    #   Target: test_gelim                                                |
    #----------------------------------------------------------------------

       test_gelim.NAME     = Test of MP File Input
       test_gelim.TYPE     = C++
       test_gelim.OBJS     = test_gelim.o
       test_gelim.DEP      = libgelim.a
       test_gelim.LDLIBS   = $(LIB_DATA_CLEANING)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

