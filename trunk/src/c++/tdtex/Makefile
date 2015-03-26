# Makefile for TDTEX
#
#  Author: Kevin Jacobs (jacobs@theopalgroup.com)
#
#  History:   1.0   kbj  Initial implementation        June  3 2002
#
#  Copyright (c) 2002  R.C. Elston
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

  TARGET_NAME = TDTEX
  TARGET      =
  TARGETS     = libtdtex.a tdtex$(EXE)
  TESTTARGETS = libtdtex.a tdtex$(EXE)
  TARPREFIX   = TDTEX
  TESTS       = runall tdtex

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS        = Configuration.cpp             \
                Parser.cpp                    \
                Transmission.cpp              \
                TransmissionTable.cpp         \
                AlleleTransmissionTable.cpp   \
                GenotypeTransmissionTable.cpp \
                ScoreCalculator.cpp           \
                Sampler.cpp                   \
                Tests.cpp

  DEP_SRCS    = tdtex.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libtdtex.a                                                |
    #----------------------------------------------------------------------

       libtdtex.a.NAME     = "Transmission-disequilibrium test library"
       libtdtex.a.TYPE     = LIB
       libtdtex.a.CXXFLAGS = -I../mped
       libtdtex.a.OBJS     = $(OBJS)
       libtdtex.a.CP       = ../lib/libtdtex.a

    #======================================================================
    #   Target: tdtex                                                     |
    #----------------------------------------------------------------------

       tdtex$(EXE).NAME     = "Transmission-disequilibrium test application"
       tdtex$(EXE).INSTALL  = yes
       tdtex$(EXE).TYPE     = C++
       tdtex$(EXE).OBJS     = tdtex.o
       tdtex$(EXE).DEP      = libtdtex.a
       tdtex$(EXE).LDLIBS   = -ltdtex $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

