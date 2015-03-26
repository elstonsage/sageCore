# Makefile for Fortran IO Objects
#
#  Authors: Kevin Jacobs (jacobs@darwin.cwru.edu)
#           Geoff Wedig  (wedig@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Apr  5 1996
#             1.01  kbj  Refinements and confusion     Jun 10 1996
#             1.1   kbj  Finally got it right          Jul  1 1996
#
#  Copyright (c) 1996  R.C. Elston

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

  TARGET_NAME = "Fortran I/O Objects"
  TARGET      =
  TARGETS     = libfortran.a 
  TESTTARGETS = libfortran.a
  VERSION     = 1.0
  TARPREFIX   = FIO
  TESTS       = true

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS      = Formatter.h   Tokenizer.h   Token_func.h

  SRCS         = Formatter.cpp Tokenizer.cpp Token_func.cpp

  OBJS         = $(SRCS:.cpp=.o)

    #======================================================================
    #   Target: libfortran.a blob                                         |
    #----------------------------------------------------------------------

       libfortran.a.NAME     = "Fortran Parsing blob"
       libfortran.a.TYPE     = LIB
       libfortran.a.OBJS     = $(OBJS)
       libfortran.a.CP       = ../lib/libfortran.a
       libfortran.a.CXXFLAGS =

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)

