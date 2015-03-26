# Makefile for IBD Generator
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

  LSFROOT = ..
  SAGEROOT= ../..

  include  $(SAGEROOT)/config/Global.make

  TARGET_NAME = IBD Generator
  TARGET      = 
  TARGETS     = libgenibd.a genibd$(EXE)
  TESTTARGETS = libgenibd.a genibd$(EXE)
  VERSION     = 1.0
  TARPREFIX   = GENIBD
  TESTS       = runall genibd

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = input.h                parser.h        params.h  \
                pair_ibd_analysis.h    peeler.h        analysis.h \
                sim_pair.h             sim_storage_ibd.h \
                sim_simulator.h        sim_ibd_analysis.h \
                genibd.h

  SRCS        = input.cpp              parser.cpp      params.cpp \
                pair_ibd_analysis.cpp  peeler.cpp      analysis.cpp \
                sim_simulator.cpp      sim_pair.cpp \
                sim_ibd_analysis.cpp

  DEP_SRCS    = genibd.cpp

  OBJS        = ${SRCS:.cpp=.o} 

    #======================================================================
    #   Target: libgenibd.a                                              |
    #----------------------------------------------------------------------

       libgenibd.a.NAME     = "IBD generator program Library"
       libgenibd.a.TYPE     = LIB
       libgenibd.a.CXXFLAGS = #-DDEBUG_VERBOSE -I. -g 
       libgenibd.a.LDFLAGS  = # --backend -non_shared  #-DDEBUG_VERBOSE -I. -g 
       libgenibd.a.OBJS     = ${OBJS}
       libgenibd.a.CP       = ../lib/libgenibd.a

    #======================================================================
    #   Target: genibd                                                   |
    #----------------------------------------------------------------------

       genibd$(EXE).NAME     = "IBD generator program"
       genibd$(EXE).INSTALL  = yes
       genibd$(EXE).TYPE     = C++
       genibd$(EXE).CXXFLAGS = #-DDEBUG_VERBOSE -I. -g 
       genibd$(EXE).LDFLAGS  = # --backend -non_shared  #-DDEBUG_VERBOSE -I. -g 
       genibd$(EXE).OBJS     = genibd.o 
       genibd$(EXE).LDLIBS   = -lgenibd $(LIB_ALL)

include  $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

