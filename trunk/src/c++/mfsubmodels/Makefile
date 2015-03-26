# Makefile for Mfsubmodels
#
#  Authors: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Oct 22 1998
#
#  Copyright (c) 1998  R.C. Elston
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

  TARGET_NAME = "Mfsubmodels"
  TARGET      =
  TARGETS     = libmfsubmodels.a
  TESTTARGETS = libmfsubmodels.a test_type_specific_submodel$(EXE) test_transformation_submodel$(EXE)

  VERSION     = 1.0
  TARPREFIX   = MFSUBMDOELS
  TESTS       = runall mfsubmodels

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  TRANSFORMATION = TransformationSubmodel.cpp TransformationParser.cpp Transformation.cpp
  TYPESPECIFIC   = TypeSubmodel.cpp TypeParser.cpp
  FPMM           = # FPMMSubmodel.cpp

  SRCS           = ${TRANSFORMATION} ${TYPESPECIFIC} ${FPMM}

  DEP_SRCS       = test_type_specific_submodel.cpp test_transformation_submodel.cpp

  OBJS           = ${SRCS:.cpp=.o}

  HEADERS        = ${SRCS:.cpp=.h}

    #======================================================================
    #   Target: libmfsubmodels.a                                               |
    #----------------------------------------------------------------------

       libmfsubmodels.a.NAME     = Maxfun Submodels
       libmfsubmodels.a.TYPE     = LIB
       libmfsubmodels.a.OBJS     = ${OBJS}
       libmfsubmodels.a.CP       = ../lib/libmfsubmodels.a

    #======================================================================
    #   Target: test_type_specific_submodel.cpp                                                     |
    #----------------------------------------------------------------------

       test_type_specific_submodel$(EXE).NAME      = Test for TypeSpecificSubmodel
       test_type_specific_submodel$(EXE).TYPE      = C++
       test_type_specific_submodel$(EXE).OBJS      = test_type_specific_submodel.o
       test_type_specific_submodel$(EXE).DEP       = libmfsubmodels.a
       test_type_specific_submodel$(EXE).LDLIBS    = $(LIB_ALL) $(LAPACK) -lm

    #======================================================================
    #   Target: test_transformation_submodel
    #----------------------------------------------------------------------

       test_transformation_submodel$(EXE).NAME      = Test for TransformationSubmodel
       test_transformation_submodel$(EXE).TYPE      = C++
       test_transformation_submodel$(EXE).OBJS      = test_transformation_submodel.o
       test_transformation_submodel$(EXE).DEP       = libmfsubmodels.a
       test_transformation_submodel$(EXE).LDLIBS    = $(LIB_ALL) $(LAPACK) -lm

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

