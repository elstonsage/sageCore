# Makefile for Filtered Pedigrees
#
#  Copyright (c) 2004  R.C. Elston
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

  TARGET_NAME = "Filtered Pedigrees and related objects" 
  TARGET      =
  TARGETS     = libfped.a
  TESTTARGETS = libfped.a filter_test$(EXE)
  TARPREFIX   = MPS
  TESTS       = runall fped


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = fped_obj.h fped_filter.cpp fped_func.h

  SRCS        = fped_obj.cpp fped_filter.cpp fped_func.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libfped.a                                             |
    #----------------------------------------------------------------------

       libfped.a.NAME     = FILTERED Pedigree Library
       libfped.a.TYPE     = LIB
       libfped.a.OBJS     = $(OBJS)
       libfped.a.CP       = ../lib/libfped.a
       libfped.a.CXXFLAGS = 

    #======================================================================
    #   Target: filter_test$(EXE)                                                 |
    #----------------------------------------------------------------------

       filter_test$(EXE).NAME     = "Multi-pedigree loop test program"
       filter_test$(EXE).TYPE     = C++
       filter_test$(EXE).OBJS     = filter_test.o
       filter_test$(EXE).DEP      = libfped.a
       filter_test$(EXE).LDLIBS   = -lfped $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

filter_test.o: ../include/fped/fped_func.ipp ../include/fped/fped_func.h


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)

