# Makefile for Numeric Routines
#
#  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial version               03-07-99
#
#  Copyright (c) 1999  R.C. Elston
#

#==========================================================================
# General target options                                                  |
#--------------------------------------------------------------------------
  LSFROOT=..
  SAGEROOT=../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = Error handler streams
  TARGET      =
  TARGETS     = libnumerics.a
  TESTTARGETS = libnumerics.a
#  TESTTARGETS = libnumerics.a numerictest
  VERSION     = 1.0
  TESTS       = true

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = cephes.h constants.h functions.h histogram.h kahan.h sinfo.h \
                mt.h corinfo.h geometric.h 

  STAT        = corinfo.cpp 

  DISTRIB     = bdtr.cpp btdtr.cpp chdtr.cpp choose.cpp fdtr.cpp gamma.cpp \
                gdtr.cpp igam.cpp igami.cpp incbet.cpp incbi.cpp nbdtr.cpp \
                ndtr.cpp ndtri.cpp pdtr.cpp stdtr.cpp unity.cpp

  LINALG      = fmatrix.cpp  # linalg.cpp 

  RNG         = mt.cpp geometric.cpp

  UTIL        = print_util.cpp

  SRCS        = ${STAT} ${DISTRIB} ${RNG} ${LINALG} $(UTIL)

  OBJS        = ${SRCS:%.cpp=%.o}
  TSRCS       = constants.cpp

  #======================================================================
  #   Target: libnumerics.a                                                |
  #----------------------------------------------------------------------

       libnumerics.a.NAME     = numerics Library
       libnumerics.a.CXXFLAGS =
       libnumerics.a.OBJS     = $(OBJS)
       libnumerics.a.SRCS     = $(SRCS)
       libnumerics.a.TYPE     = LIB
       libnumerics.a.CP       = ../lib/libnumerics.a



  #======================================================================
  #   Target: constants                                                 |
  #----------------------------------------------------------------------

       constants$(EXE).NAME      = 
       constants$(EXE).TYPE      = C++
       constants$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE
       constants$(EXE).OBJS      = constants.o 
       constants$(EXE).LDFLAGS   = -L.
       constants$(EXE).LDLIBS    = $(LIB_ALL)

  #======================================================================
  #   Target: testcor                                                 |
  #----------------------------------------------------------------------

       testcor$(EXE).NAME      = 
       testcor$(EXE).TYPE      = C++
       testcor$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE
       testcor$(EXE).OBJS      = testcor.o 
       testcor$(EXE).LDFLAGS   = -L.
       testcor$(EXE).LDLIBS    = $(LIB_ALL)

  #======================================================================
  #   Target: test_log_double                                                 |
  #----------------------------------------------------------------------

       test_log_double$(EXE).NAME      = 
       test_log_double$(EXE).TYPE      = C++
       test_log_double$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE
       test_log_double$(EXE).OBJS      = test_log_double.o 
       test_log_double$(EXE).LDFLAGS   = -L.
       test_log_double$(EXE).LDLIBS    = $(LIB_ALL)

  #======================================================================
  #   Target: sinfo                                                     |
  #----------------------------------------------------------------------

       sinfo$(EXE).NAME      = 
       sinfo$(EXE).TYPE      = C++
       sinfo$(EXE).CXXFLAGS  = 
       sinfo$(EXE).OBJS      = sinfo.o 
       sinfo$(EXE).LDFLAGS   = -L.
       sinfo$(EXE).LDLIBS    = $(LIB_ALL)

  #======================================================================
  #   Target: svdtest                                                   |
  #----------------------------------------------------------------------

       svdtest$(EXE).NAME      = 
       svdtest$(EXE).TYPE      = C++
       svdtest$(EXE).CXXFLAGS  = 
       svdtest$(EXE).OBJS      = svdtest.o 
       svdtest$(EXE).LDFLAGS   = -L.
       svdtest$(EXE).LDLIBS    = $(LIB_ALL)

  #======================================================================
  #   Target: fmatrixtest                                               |
  #----------------------------------------------------------------------

       fmatrixtest$(EXE).NAME      = 
       fmatrixtest$(EXE).TYPE      = C++
       fmatrixtest$(EXE).CXXFLAGS  = 
       fmatrixtest$(EXE).OBJS      = fmatrixtest.o 
       fmatrixtest$(EXE).LDFLAGS   = -L.
       fmatrixtest$(EXE).LDLIBS    = $(LIB_ALL)

  #======================================================================
  #   Target: linalgtest                                                |
  #----------------------------------------------------------------------

       linalgtest$(EXE).NAME      = 
       linalgtest$(EXE).TYPE      = C++
       linalgtest$(EXE).CXXFLAGS  = 
       linalgtest$(EXE).OBJS      = linalgtest.o 
       linalgtest$(EXE).LDFLAGS   = -L.
       linalgtest$(EXE).LDLIBS    = $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)

