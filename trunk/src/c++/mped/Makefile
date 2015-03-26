# Makefile for New Multipedigree
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

  TARGET_NAME = New MultiPedigree
  TARGET      =
  TARGETS     = libmped.a 
  TESTTARGETS = libmped.a a2_test$(EXE) # mpfiletest
  VERSION     = 1.0
  TESTS       = runall mped


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS      = mp.h mpiterator.h  spbaseiter.h  spiterator.h   \
                 mpbase.h sp.h sptypes.h mpfwd.h                 \
                 spbase.h spiter.h

  SRCS         = \
                 spbase0.cpp          \
                 member_base.cpp      \
                 family_base.cpp      \
                 subpedigree_base.cpp \
                 pedigree_base.cpp    \
                 PedigreeBuilder.cpp  \
                 spbaseiter.cpp       \
                 mpbase.cpp           \
                 mp_utilities.cpp

  DEP_SRCS     = a2_test.cpp

  OBJS         = $(SRCS:.cpp=.o)

    #======================================================================
    #   Target: libmped.a                                                 |
    #----------------------------------------------------------------------

       libmped.a.NAME     = New MultiPedigree Library
       libmped.a.TYPE     = LIB
       libmped.a.OBJS     = $(OBJS)
       libmped.a.CP       = ../lib/libmped.a
       libmped.a.CXXFLAGS = -I.

    #======================================================================
    #   Target: a2_test                                                   |
    #----------------------------------------------------------------------

       a2_test$(EXE).NAME     = New MultiPedigree Test
       a2_test$(EXE).TYPE     = C++
       a2_test$(EXE).OBJS     = a2_test.o
       a2_test$(EXE).LDLIBS   = ../lib/libmped.a
       a2_test$(EXE).DEP      = libmped.a
       a2_test$(EXE).CXXFLAGS = -L.


include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

