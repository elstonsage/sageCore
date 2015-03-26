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

  TARGET_NAME = "Decipher"
  TARGET      = 
  TARGETS     = libdecipher.a decipher$(EXE)  
  TESTTARGETS = libdecipher.a decipher$(EXE) 
  TARPREFIX   = DECIPHER
  TESTS       = runall decipher

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = input.h em.h instructions.h \
                parser.h analysis.h unrelated.h \
                family.h output.h rebuilt.h \
                sub_pop_shuffler.h founders.h pool.h \
                definitions.h partitioning.h

  SRCS        = input.cpp  em.cpp instructions.cpp \
                parser.cpp analysis.cpp unrelated.cpp \
                family.cpp output.cpp rebuilt.cpp \
                sub_pop_shuffler.cpp founders.cpp pool.cpp \
                definitions.cpp partitioning.cpp
  
  DEP_SRCS    = decipher.cpp

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS = ${SRCS:.cpp=.h,}

    #======================================================================
    #   Target: libdecipher.a                                             |
    #----------------------------------------------------------------------

       libdecipher.a.NAME     = "decipher library"
       libdecipher.a.TYPE     = LIB
       libdecipher.a.CXXFLAGS = -I../mped 
       libdecipher.a.OBJS     = $(OBJS)
       libdecipher.a.CP       = ../lib/libdecipher.a


    #======================================================================
    #   Target: decipher                                                  |
    #----------------------------------------------------------------------

       decipher$(EXE).NAME     = DECIPHER Application
       decipher$(EXE).INSTALL  = yes
       decipher$(EXE).TYPE     = C++
       decipher$(EXE).OBJS     = decipher.o 
       decipher$(EXE).DEP      = libdecipher.a
       decipher$(EXE).LDLIBS   = -ldecipher $(LIB_ALL)
                         
  		                        
    #======================================================================
    #   Target: scratch                                                   |
    #----------------------------------------------------------------------

       scratch$(EXE).NAME     = scratch 
       scratch$(EXE).TYPE     = C++
       scratch$(EXE).OBJS     = scratch.o
       scratch$(EXE).LDLIBS   = -ldecipher  $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

