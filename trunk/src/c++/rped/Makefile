# Makefile for New General Referenced Pedigrees
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

  TARGET_NAME = "New General Referenced Pedigrees and related objects" 
  TARGET      =
  TARGETS     = librped.a
  TESTTARGETS = librped.a test_rp_info mpfiletest loop_test
  TARPREFIX   = MPS
  TESTS       = runall rped


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = \
                genome_description.h \
                Invalidator.h        \
                loop.h               \
                rped.h               \
                rpfile.h             \

  SRCS        = \
                genome_description.cpp \
                loop.cpp               \
                rped.cpp               \
                rpfile.cpp             \
                rpfile_fortran.cpp     \
                rpfile_delimited.cpp   \

  DEP_SRCS    = test_rp_info.cpp mpfiletest.cpp loop_test.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: librped.a                                                 |
    #----------------------------------------------------------------------

       librped.a.NAME     = New General Referenced Pedigree Library
       librped.a.TYPE     = LIB
       librped.a.OBJS     = $(OBJS)
       librped.a.CP       = ../lib/librped.a
       librped.a.CXXFLAGS = 


    #======================================================================
    #   Target: loop_test                                                 |
    #----------------------------------------------------------------------

       loop_test.NAME     = "Multi-pedigree loop test program"
       loop_test.TYPE     = C++
       loop_test.DEP      = librped.a
       loop_test.OBJS     = loop_test.o
       loop_test.LDLIBS   = $(LIB_PEDIGREE_DATA)

    #======================================================================
    #   Target: mpfiletest                                                |
    #----------------------------------------------------------------------

       mpfiletest.NAME     = Test of MP File Input
       mpfiletest.TYPE     = C++
       mpfiletest.OBJS     = mpfiletest.o 
       mpfiletest.DEP      = librped.a
       mpfiletest.LDLIBS   = $(LIB_PEDIGREE_DATA)

    #======================================================================
    #   Target: test_rp_info                                              |
    #----------------------------------------------------------------------

       test_rp_info.NAME     = Test of MP File Input
       test_rp_info.TYPE     = C++
       test_rp_info.OBJS     = test_rp_info.o 
       test_rp_info.DEP      = librped.a
       test_rp_info.LDLIBS   = $(LIB_PEDIGREE_DATA)


include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

