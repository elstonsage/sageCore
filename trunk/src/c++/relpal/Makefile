# Makefile for Relpal
#
#  Authors: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Oct 22 1998
#
#  Copyright (c) 1998  R.C. Elston
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

  LSFROOT = ..
  SAGEROOT= ../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = Relpal
  TARGET      =
  TARGETS     = librelpal.a relpal$(EXE)
  TESTTARGETS = librelpal.a relpal$(EXE)
  VERSION     = 1.0
  TARPREFIX   = RELPAL
  TESTS       = runall relpal

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = definitions.h                                       \
                input.h                   params.h                  \
                model.h                   parser.h                  \
                semi_definite.h           covariance_calculator.h   \
                two_level_calculator.h    two_level_test.h          \
                analysis.h                output.h                  \
                relpal.h

  SRCS        = input.cpp                 params.cpp                \
                model.cpp                 parser.cpp                \
                semi_definite.cpp         covariance_calculator.cpp \
                calc_base.cpp             calc_gls.cpp              \
                calc_score.cpp            two_level_test_base.cpp   \
                two_level_regression.cpp  two_level_score_test.cpp  \
                analysis.cpp              output.cpp

  DEP_SRCS    = relpal.cpp                test_matrix.cpp  test_pibd.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: librelpal.a                                               |
    #----------------------------------------------------------------------

       librelpal.a.NAME     = "Relpal Library"
       librelpal.a.TYPE     = LIB
       librelpal.a.OBJS     = $(OBJS)
       librelpal.a.CP       = ../lib/librelpal.a

    #======================================================================
    #   Target: relpal                                                    |
    #----------------------------------------------------------------------

       relpal$(EXE).NAME      = Relpal
       relpal$(EXE).INSTALL   = yes
       relpal$(EXE).TYPE      = C++
       relpal$(EXE).OBJS      = relpal.o
       relpal$(EXE).DEP       = librelpal.a
       relpal$(EXE).LDLIBS    = -lrelpal -lpalbase $(LIB_ALL) 

    #======================================================================
    #   Target: test_matrix                                               |
    #----------------------------------------------------------------------

       test_matrix$(EXE).NAME      = Relpal
       test_matrix$(EXE).INSTALL   = yes
       test_matrix$(EXE).TYPE      = C++
       test_matrix$(EXE).OBJS      = test_matrix.o
       test_matrix$(EXE).DEP       = librelpal.a
       test_matrix$(EXE).LDLIBS    = -lpalbase $(LIB_ALL) 

    #======================================================================
    #   Target: test_pibd                                                 |
    #----------------------------------------------------------------------

       test_pibd$(EXE).NAME      = "test prior ibd"
       test_pibd$(EXE).INSTALL   = yes
       test_pibd$(EXE).TYPE      = C++
       test_pibd$(EXE).OBJS      = input.o test_pibd.o
       test_pibd$(EXE).DEP       = librelpal.a
       test_pibd$(EXE).LDLIBS    = -lpalbase $(LIB_ALL) 

include $(SAGEROOT)/config/Rules.make



# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

