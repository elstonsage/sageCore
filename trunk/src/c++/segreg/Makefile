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

  TARGET_NAME = "Regression analysis"
  TARGET      =
  TARGETS     = libsegreg.a segreg$(EXE)
  TESTTARGETS = libsegreg.a segreg$(EXE) test_type_description$(EXE) test_segreg$(EXE) test_parser$(EXE)
  TARPREFIX   = SEGREG
  TESTS       = runall segreg

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  NEW         = types/TypeDescription.cpp

  MODELS      = type_sub_model.cpp                                         \
                fpmm_sub_model.cpp freq_sub_model.cpp                      \
                transmission_sub_model.cpp onset_sub_model.cpp             \
                ascertainment_sub_model.cpp cov_sub_model.cpp              \
                resid_sub_model.cpp prev_sub_model.cpp model.cpp

  DATA_MGRS   = PedigreeDataSet.cpp

  CALC        = member_calculator.cpp SL_calculator.cpp                    \
                regressive_peeler.cpp FPMM_peeler.cpp mlm_peeler.cpp       \
                MlmCorrelationVerifier.cpp                                 \
                segreg_calculator.cpp mean_split.cpp                       \
                freq_estimates.cpp polygenic_transition_calculator.cpp     \
                iterative_models.cpp mlm_resid_corr.cpp                    \
                RegPenetranceCommonAdjustments.cpp RegPenetranceContext.cpp\
                RegPenetranceCalculator.cpp

  CONTROL     = parser.cpp segreg_input.cpp segreg_utilities.cpp           \
		pen_function.cpp analysis.cpp AnalysisController.cpp       \
                analysis_out.cpp

  SRCS        = $(NEW) $(MODELS) $(DATA_MGRS) $(CALC) $(CONTROL)

  DEP_SRCS   = segreg.cpp test_segreg.cpp test_parser.cpp types/test_type_description.cpp

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS = ${SRCS:.cpp=.h,}

    #======================================================================
    #   Target: libsegreg.a                                              |
    #----------------------------------------------------------------------

       libsegreg.a.NAME     = "segreg library"
       libsegreg.a.TYPE     = LIB
       libsegreg.a.CXXFLAGS = -I../mped
       libsegreg.a.SRCS     = $(SRCS)
       libsegreg.a.OBJS     = $(OBJS)
       libsegreg.a.CP       = ../lib/libsegreg.a

    #======================================================================
    #   Target: segreg
    #----------------------------------------------------------------------

       segreg$(EXE).NAME     = SEGREG Application
       segreg$(EXE).INSTALL  = yes
       segreg$(EXE).TYPE     = C++
       segreg$(EXE).OBJS     = segreg.o
       segreg$(EXE).DEP      = libsegreg.a
       segreg$(EXE).LDLIBS   = -lsegreg  $(LIB_ALL) 

    #======================================================================
    #   Target: test_type_description
    #----------------------------------------------------------------------

       test_type_description$(EXE).NAME     = SEGREG Application
       test_type_description$(EXE).INSTALL  = yes
       test_type_description$(EXE).TYPE     = C++
       test_type_description$(EXE).OBJS     = types/test_type_description.o
       test_type_description$(EXE).DEP      =
       test_type_description$(EXE).LDLIBS   = -lsegreg  $(LIB_ALL) 
    #======================================================================
    #   Target: testof
    #----------------------------------------------------------------------

       testof.NAME     = SEGREG Application
       testof.TYPE     = C++
       testof.OBJS     = testof.o

    #======================================================================
    #   Target: test_parser                                                |
    #----------------------------------------------------------------------

       test_parser$(EXE).NAME     = Test of segreg parser
       test_parser$(EXE).TYPE     = C++
       test_parser$(EXE).OBJS     = test_parser.o
       test_parser$(EXE).DEP      = libsegreg.a
       test_parser$(EXE).LDLIBS   = -lsegreg  $(LIB_ALL) 


    #======================================================================
    #   Target: test_segreg
    #----------------------------------------------------------------------

       test_segreg$(EXE).NAME     = Test of segreg calculators
       test_segreg$(EXE).TYPE     = C++
       test_segreg$(EXE).OBJS     = test_segreg.o
       test_segreg$(EXE).DEP      = libsegreg.a
       test_segreg$(EXE).LDLIBS   = -lsegreg  $(LIB_ALL) 

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

