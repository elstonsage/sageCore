# Makefile for Error Handler
#
#  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial version               03-09-98
#
#  Copyright (c) 1998  R.C. Elston
#

#==========================================================================
# General target options                                                  |
#--------------------------------------------------------------------------
  LSFROOT=..
  SAGEROOT=../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = Error handler streams
  TARGET      =
  TARGETS     = liberror.a #liberror_reg.a
  TESTTARGETS = liberror.a # errproc$(EXE) # errortest$(EXE)
  VERSION     = 1.0
  TESTS       = true

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = errorbuf.h errormanip.h errorstream.h \
                bufferederrorbuf.h bufferederrorstream.h
  SRCS        = errorstream.cpp
  OBJS        = ${SRCS:%.cpp=%.o}
  TSRCS       = errortest.cpp

  #======================================================================
  #   Target: liberror.a                                                |
  #----------------------------------------------------------------------

       liberror.a.NAME     = error Library
       liberror.a.OBJS     = $(OBJS)
       liberror.a.SRCS     = $(SRCS)
       liberror.a.TYPE     = LIB
       liberror.a.CP       = ../lib/liberror.a

  #======================================================================
  #   Target: liberror_reg.a                                                |
  #----------------------------------------------------------------------

       liberror_reg.a.NAME     = error registry Library
       liberror_reg.a.OBJS     = Registry.o RegistryRenderer.o
       liberror_reg.a.SRCS     = Registry.cpp RegistryRenderer.cpp
       liberror_reg.a.TYPE     = LIB
       liberror_reg.a.CXXFLAGS = -O0
       liberror_reg.a.CP       = ../lib/liberror_reg.a

  #======================================================================
  #   Target: errortest$(EXE)                                                 |
  #----------------------------------------------------------------------

       errortest$(EXE).NAME      = Error Handler Test Program
       errortest$(EXE).TYPE      = C++
       errortest$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE 
       errortest$(EXE).OBJS      = errortest.o 
       errortest$(EXE).LDFLAGS   = -L.
       errortest$(EXE).LDLIBS    = -lerror $(LIB_EXT_APIS)

  #======================================================================
  #   Target: registry_test$(EXE)                                                 |
  #----------------------------------------------------------------------

       registry_test$(EXE).NAME      = ErrorRegistry Test Program
       registry_test$(EXE).TYPE      = C++
       registry_test$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE
       registry_test$(EXE).OBJS      = registry_test.o 
       registry_test$(EXE).LDFLAGS   = -L.
       registry_test$(EXE).DEP       = liberror.a
       registry_test$(EXE).LDLIBS    = $(LIB_ALL)

  #======================================================================
  #   Target: errproc$(EXE)                                                 |
  #----------------------------------------------------------------------

       errproc$(EXE).NAME      = ErrorRegistry-to-Latex converter program
       errproc$(EXE).TYPE      = C++
       errproc$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE
       errproc$(EXE).OBJS      = errproc.o 
       errproc$(EXE).LDFLAGS   = -L.
       errproc$(EXE).DEP       = liberror.a
       errproc$(EXE).LDLIBS    = $(LIB_ALL)

  #======================================================================
  #   Target: buffertest                                                 |
  #----------------------------------------------------------------------

       buffertest$(EXE).NAME      = Error Handler Test Program
       buffertest$(EXE).TYPE      = C++
       buffertest$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE
       buffertest$(EXE).OBJS      = buffertest.o 
       buffertest$(EXE).LDFLAGS   = -L.
       buffertest$(EXE).LDLIBS    = $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)

