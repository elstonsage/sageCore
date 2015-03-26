# Makefile for function objects.
#
#  Author: Dan Baechle
#
#  History:   Created 12/11/00
#
#  Copyright (c) 2000  R.C. Elston
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

  TARGET_NAME = MultiPedigree variable transformation
  TARGET      =
  TARGETS     = libfunc.a 
  TESTTARGETS = libfunc.a test_func test_parse test_expression
  TARPREFIX   = FUNC
  TESTS       = runall func


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS        = Function.cpp PythonInterface.cpp FunctionParser.cpp evalfunc.cpp \
                AdjustedTraitCreator.cpp tai.cpp twp.cpp Expression.cpp \
                ParentOfOrigin.cpp

  DEP_SRCS    = test_func.cpp test_parse.cpp test_expression.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libfunc.a                                                |
    #----------------------------------------------------------------------

       libfunc.a.NAME     = Function library
       libfunc.a.CXXFLAGS =
       libfunc.a.TYPE     = LIB
       libfunc.a.OBJS     = $(OBJS)
       libfunc.a.CP       = ../lib/libfunc.a


    #======================================================================
    #   Target: test_func                                                 |
    #----------------------------------------------------------------------

       test_func.NAME     = Test of function code
       test_func.TYPE     = C++
       test_func.OBJS     = test_func.o 
       test_func.DEP      = libfunc.a
       test_func.LDLIBS   = $(LIB_DATA_CLEANING)

    #======================================================================
    #   Target: test_expression                                                 |
    #----------------------------------------------------------------------

       test_expression.NAME     = Test of expression code
       test_expression.TYPE     = C++
       test_expression.OBJS     = test_expression.o 
       test_expression.DEP      = libfunc.a
       test_expression.LDLIBS   = $(LIB_DATA_CLEANING)

    #======================================================================
    #   Target: test_parse                                                |
    #----------------------------------------------------------------------

       test_parse.NAME     = Test of parsing code
       test_parse.TYPE     = C++
       test_parse.OBJS     = test_parse.o 
       test_parse.DEP      = libfunc.a
       test_parse.LDLIBS   = $(LIB_DATA_CLEANING)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

