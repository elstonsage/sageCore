# Makefile for New Multi Locus
#
#  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial version            2 27 1998
#
#  Copyright (c) 1998  R.C. Elston

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

  TARGET_NAME = New Multi Locus
  TARGET      =
  TARGETS     = libmlocus.a 
  TESTTARGETS = libmlocus.a test_mlocus$(EXE) 
  VERSION     = 1.0
  TESTS       = runall mlocus

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS = genotype.cpp phenotype.cpp phmodel.cpp penmodel.cpp imodel.cpp mfile.cpp

  DEP_SRCS = test_mlocus.cpp test_phmodel.cpp test_penmodel.cpp test_imodel.cpp test_mfile.cpp

  OBJS = $(SRCS:.cpp=.o)

    #======================================================================
    #   Target: libmlocus.a blob                                          |
    #----------------------------------------------------------------------

       libmlocus.a.NAME     = New Multi Locus Library
       libmlocus.a.TYPE     = LIB
       libmlocus.a.OBJS     = $(OBJS)
       libmlocus.a.CP       = ../lib/libmlocus.a
       libmlocus.a.CXXFLAGS = -I.

    #======================================================================
    #   Target: test_mlocus                                               |
    #----------------------------------------------------------------------

       test_mlocus$(EXE).NAME     = New Multi Locus Test
       test_mlocus$(EXE).TYPE     = C++
       test_mlocus$(EXE).DEP      = libmlocus.a
       test_mlocus$(EXE).OBJS     = test_mlocus.o test_phmodel.o test_penmodel.o test_imodel.o test_mfile.o
       test_mlocus$(EXE).LDLIBS   = -lmlocus -loutput $(LIB_CORE)
       test_mlocus$(EXE).CXXFLAGS = -L.

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

