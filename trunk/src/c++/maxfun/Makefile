# Makefile for Maxfun
#
#  Authors: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Oct 22 1998
#
#  Copyright (c) 1998  R.C. Elston
#
#  NOTE:  All S.A.G.E. programmers; this Makefile is a good example of how
#         to derive classes from LSF and build them.

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

  TARGET_NAME = "Maxfun"
  TARGET      =
  TARGETS     = libmaxfun.a
  TESTTARGETS = libmaxfun.a maxex2$(EXE) maxex3$(EXE) maxex4$(EXE)
  TARPREFIX   = MAXFUN
  TESTS       = true

#  TESTS       = cd tests/test1; rm maxex.out 2> /dev/null;  \
#                for i in 1 2 3 4 5 6; do              \
#                  ../../maxex2$(EXE) $$i >> maxex.out;          \
#                done;                                 \
#                diff -u maxex.exp maxex.out > diffs

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS = maxfun.cpp cmaxfun.cpp sub_model.cpp transf_sub_model.cpp maxfun_utilities.cpp

  DEP_SRCS = maxex2.cpp maxex3.cpp maxex4.cpp

  OBJS = ${SRCS:%.cpp=%.o}

    #======================================================================
    #   Target: libmaxfun.a                                               |
    #----------------------------------------------------------------------

       libmaxfun.a.NAME     = Maxfun Library
       libmaxfun.a.TYPE     = LIB
       libmaxfun.a.OBJS     = ${OBJS}
       libmaxfun.a.CP       = ../lib/libmaxfun.a
       libmaxfun.a.LDLIBS   = $(LIB_ALL)

    #======================================================================
    #   Target: maxex                                                     |
    #----------------------------------------------------------------------

       maxex$(EXE).NAME      = Simple Maxfun example
       maxex$(EXE).INSTALL   = no
       maxex$(EXE).TYPE      = C++
       maxex$(EXE).SRCS      = maxex.cpp
       maxex$(EXE).OBJS      = maxex.o
       maxex$(EXE).DEP       = libmaxfun.a
       maxex$(EXE).LDLIBS    = $(LIB_ALL)

    #======================================================================
    #   Target: maxex2                                                     |
    #----------------------------------------------------------------------

       maxex2$(EXE).NAME      = Simple Maxfun example
       maxex2$(EXE).INSTALL   = no
       maxex2$(EXE).TYPE      = C++
       maxex2$(EXE).SRCS      = maxex2.cpp
       maxex2$(EXE).OBJS      = maxex2.o
       maxex2$(EXE).LDLIBS    = $(LIB_ALL)

    #======================================================================
    #   Target: maxex3                                                     |
    #----------------------------------------------------------------------

       maxex3$(EXE).NAME      = Simple Maxfun example
       maxex3$(EXE).INSTALL   = no
       maxex3$(EXE).TYPE      = C++
       maxex3$(EXE).SRCS      = maxex3.cpp
       maxex3$(EXE).OBJS      = maxex3.o
       maxex3$(EXE).LDLIBS    = $(LIB_ALL)

    #======================================================================
    #   Target: maxex4                                                     |
    #----------------------------------------------------------------------

       maxex4$(EXE).NAME      = Simple Maxfun example
       maxex4$(EXE).INSTALL   = no
       maxex4$(EXE).TYPE      = C++
       maxex4$(EXE).SRCS      = maxex4.cpp
       maxex4$(EXE).OBJS      = maxex4.o
       maxex4$(EXE).LDLIBS    = $(LIB_ALL)

    #======================================================================
    #   Target: maxex_opt                                                     |
    #----------------------------------------------------------------------

       maxex_opt.NAME      = Simple Maxfun example
       maxex_opt.INSTALL   = no
       maxex_opt.TYPE      = C++
       maxex_opt.SRCS      = maxex_opt.cpp
       maxex_opt.OBJS      = maxex_opt.o
       maxex_opt.DEP       = libmaxfunapi.a
       maxex_opt.LDLIBS    = $(LIB_ALL)

    #======================================================================
    #   Target: maxtest                                                     |
    #----------------------------------------------------------------------

       maxtest.NAME      = Maxfun API Test
       maxtest.INSTALL   = yes
       maxtest.TYPE      = C++
       maxtest.SRCS      = maxtest.cpp
       maxtest.OBJS      = maxtest.o
       maxtest.DEP       = libmaxfunapi.a
       maxtest.LDLIBS    = $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

