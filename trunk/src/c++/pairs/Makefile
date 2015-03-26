# Makefile for New Pairs and Pair Generation
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

  TARGET_NAME = "New Relative Pairs and related objects"
  TARGET      =
  TARGETS     = libpairs.a 
  TESTTARGETS = libpairs.a testrelmatrix$(EXE) testrelpair$(EXE) testrelpair_simple$(EXE) testindfilter$(EXE)
  TARPREFIX   = PAIRSNEW
  TESTS       = runall pairs


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = reltype.h   reltypename.h   relmatrix.h   stringbank.h relpair.h

  SRCS        = reltype.cpp reltypename.cpp relmatrix.cpp stringbank.cpp \
                long_reltypename.cpp relpair.cpp

  DEP_SRCS    = testrelmatrix.cpp testrelpair.cpp testrelpair_simple.cpp \
                testindfilter.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libpairs.a                                                |
    #----------------------------------------------------------------------

       libpairs.a.NAME     = New Pairs Library
       libpairs.a.TYPE     = LIB
       libpairs.a.OBJS     = $(OBJS)
       libpairs.a.CP       = ../lib/libpairs.a


    #======================================================================
    #   Target: testrelmatrix                                             |
    #----------------------------------------------------------------------

       testrelmatrix$(EXE).NAME     = Test of relationship matrix
       testrelmatrix$(EXE).TYPE     = C++
       testrelmatrix$(EXE).OBJS     = testrelmatrix.o 

	# Please note: Flexlm is double-listed at the end, because of very very very bizarre
	# linker problems. Specifically, when it is *not* double-listed, the linker reports a bunch
	# of unresolved symbols from the multipedigree. I'm still working on coming up with an elegant
	# solution.

       testrelmatrix$(EXE).LDLIBS   = $(LIB_PEDIGREE_DATA)
       testrelmatrix$(EXE).DEP      = libpairs.a

    #======================================================================
    #   Target: testrelpair_simple                                        |
    #----------------------------------------------------------------------

       testrelpair_simple$(EXE).NAME     = Simple test of pair generator
       testrelpair_simple$(EXE).TYPE     = C++
       testrelpair_simple$(EXE).OBJS     = testrelpair_simple.o 
       testrelpair_simple$(EXE).LDLIBS   = $(LIB_PEDIGREE_DATA)
       testrelpair_simple$(EXE).DEP      = libpairs.a
       
       
    #======================================================================
    #   Target: testrelpair                                               |
    #----------------------------------------------------------------------

       testrelpair$(EXE).NAME     = Test of filtering pair generator
       testrelpair$(EXE).TYPE     = C++
       testrelpair$(EXE).OBJS     = testrelpair.o 
       testrelpair$(EXE).LDLIBS   = $(LIB_PEDIGREE_DATA)
       testrelpair$(EXE).DEP      = libpairs.a
       
    #======================================================================
    #   Target: testindfilter                                             |
    #----------------------------------------------------------------------

       testindfilter$(EXE).NAME     = Test of individual filter
       testindfilter$(EXE).TYPE     = C++
       testindfilter$(EXE).OBJS     = testindfilter.o 
       testindfilter$(EXE).LDLIBS   = $(LIB_PEDIGREE_DATA)
       testindfilter$(EXE).DEP      = libpairs.a


include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

