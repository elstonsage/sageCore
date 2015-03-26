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
  TARGETS     = libmlod.a mlod$(EXE)
  TESTTARGETS = libmlod.a test_parser$(EXE) mlod$(EXE)
  TARPREFIX   = MLOD.NEW
  TESTS       = runall mlod

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  HEADERS     = analysis_parameters.h    \
                analysis.h               \
                parser.h                 \
                data.h                   \
                PedigreeAnalysisSample.h \
                lod_table.h

  SRCS        = analysis_parameters.cpp    \
                parser.cpp                 \
                data.cpp                   \
                PedigreeAnalysisSample.cpp \
                lod_table.cpp              \
                analysis_results.cpp       \
                analysis_data.cpp          \
                lod_score_analyzer.cpp     \
                analysis_out.cpp           \
                analysis.cpp

  DEP_SRCS    = mlod.cpp test_parser.cpp

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS = ${SRCS:.cpp=.h,}

    #======================================================================
    #   Target: libmlod.a
    #----------------------------------------------------------------------

       libmlod.a.NAME     = "mlod library"
       libmlod.a.TYPE     = LIB
       libmlod.a.CXXFLAGS = -I../mped
       libmlod.a.OBJS     = $(OBJS)
       libmlod.a.CP       = ../lib/libmlod.a

    #======================================================================
    #   Target: test_parser                                                |
    #----------------------------------------------------------------------

       test_parser$(EXE).NAME     = Test of mlod parser
       test_parser$(EXE).TYPE     = C++
       test_parser$(EXE).OBJS     = test_parser.o
       test_parser$(EXE).DEP      = libmlod.a
       test_parser$(EXE).LDLIBS   = -lmlod $(LIB_ALL) 

    #======================================================================
    #   Target: mlod                                                |
    #----------------------------------------------------------------------

       mlod$(EXE).NAME     = Test of mlod parser
       mlod$(EXE).TYPE     = C++
       mlod$(EXE).OBJS     = mlod.o
       mlod$(EXE).DEP      = libmlod.a
       mlod$(EXE).LDLIBS   = -lmlod $(LIB_ALL)
       mlod$(EXE).INSTALL  = yes

       example.NAME     = Test of mlod parser
       example.TYPE     = C++
       example.OBJS     = math.o


include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

