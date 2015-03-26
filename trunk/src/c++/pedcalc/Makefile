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
# General target options                                                  
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

  TARGET_NAME = 
  TARGET      = 
  TARGETS     = 
  TESTTARGETS = test_bin_pen_calc$(EXE) test_fam_resid_adj$(EXE)
  TESTS       = runall ped_calc

#==========================================================================
# Source/object file lists                                                
#--------------------------------------------------------------------------

    #======================================================================
    #   Target: tests                                                |
    #----------------------------------------------------------------------

       test_fam_resid_adj$(EXE).NAME     = Test of familial_resid_adjustment objects
       test_fam_resid_adj$(EXE).TYPE     = C++
       test_fam_resid_adj$(EXE).OBJS     = test_fam_resid_adj.o 
       test_fam_resid_adj$(EXE).DEP      = 
       test_fam_resid_adj$(EXE).LDLIBS   = $(LIB_PEDIGREE_DATA) 
       test_fam_resid_adj$(EXE).INSTALL  = yes

       test_bin_pen_calc$(EXE).NAME     = Test of binary penetrance calculator object
       test_bin_pen_calc$(EXE).TYPE     = C++
       test_bin_pen_calc$(EXE).OBJS     = test_bin_pen_calc.o 
       test_bin_pen_calc$(EXE).DEP      = 
       test_bin_pen_calc$(EXE).LDLIBS   = $(LIB_TOOLS)
       test_bin_pen_calc$(EXE).INSTALL  = yes


include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

