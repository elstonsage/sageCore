# Makefile for LSF 
#
#  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Apr  5 1996
#             1.01  kbj  Refinements and confusion     Jun 10 1996
#             1.1   kbj  Finally got it right
#
#  Copyright (c) 1996  R.C. Elston
#

#==========================================================================
# General target options                                                  |
#--------------------------------------------------------------------------
  LSFROOT=..
  SAGEROOT=../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = LSF and related objects
  TARGET      =
  TARGETS     = libLSF.a
  TESTTARGETS = libLSF.a LSFtest$(EXE) xml_convert$(EXE)
  VERSION     = 1.75
  TARPREFIX   = LSF
  TESTS       = runall LSF

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HDRS        = NameMgr.h   Attr.h  LSF.h  LSFfile.h  parse_ops.h  \
                LSFtypes.h  LSFfactory.h   LSFinit.h  LSFsymbol.h  \
                ErrorHandler.h LSFiter.h LSFexpr.h var_ops.h SymbolInput.h \
                SAGEapp.h

  SRCS        = Attr.cpp ErrorHandler.cpp LSF.cpp LSFexpr.cpp        \
                LSFfactory.cpp LSFfile.cpp LSFiter.cpp LSFsymbol.cpp \
                LSFtypes.cpp NameMgr.cpp SymbolInput.cpp \
                parse_ops.cpp var_ops.cpp XMLConverter.cpp

  DEP_SRCS    = LSFtest.cpp xml_convert.cpp

  OBJS        = ${SRCS:%.cpp=%.o}

  #======================================================================
  #   Target: libLSF.a                                                  |
  #----------------------------------------------------------------------

       libLSF.a.NAME     = LSF Library
       libLSF.a.CXXFLAGS = 
       libLSF.a.OBJS     = $(OBJS)
       libLSF.a.SRCS     = $(SRCS)
       libLSF.a.TYPE     = LIB
       libLSF.a.CP       = ../lib/libLSF.a

  #======================================================================
  #   Target: LSFtest$(EXE)                                             |
  #----------------------------------------------------------------------

       LSFtest$(EXE).NAME      = LSF test program
       LSFtest$(EXE).TYPE      = C++
       LSFtest$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE
       LSFtest$(EXE).OBJS      = LSFtest.o 
       LSFtest$(EXE).LDLIBS   = libLSF.a -lerror
       LSFtest$(EXE).DEP      = libLSF.a

  #======================================================================
  #   Target: xml_convert$(EXE)                                         |
  #----------------------------------------------------------------------

       xml_convert$(EXE).NAME      = LSF-to-XML converter
       xml_convert$(EXE).TYPE      = C++
       xml_convert$(EXE).CXXFLAGS  = -DDEBUG_VERBOSE
       xml_convert$(EXE).OBJS      = xml_convert.o 
       xml_convert$(EXE).LDLIBS   = libLSF.a -lerror
       xml_convert$(EXE).DEP      = libLSF.a 

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)


