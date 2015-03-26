# Makefile for Multipoint/MCMC Objects
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

  SAGEROOT=../..
  LSFROOT=..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = "Inheritance Vectors"
  TARGETS     = liblvec.a
  TESTTARGETS = liblvec.a test_ivg test_nodes
  VERSION     = 1.0
  TARPREFIX   = IV
  TESTS       = runall lvec

#==========================================================================
# Source/object file lists                                                |

#--------------------------------------------------------------------------

  HEADERS     = meiosis_map.h   inheritance_vector.h   node.h \
                iv_generator.h  codom_ivgen.h dgraph.h lvec_allocator.h \
                fft_bit_count.h lvector.h              mpoint_like.h \
                fixed_bit_calculator.h

  SRCS        = ${HEADERS:.h=.cpp}

  DEP_SRCS    = test_ivg.cpp test_ivg_input.cpp test_nodes.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: liblvec.a                                                 |
    #----------------------------------------------------------------------

       liblvec.a.NAME     = "Likelihood Vector Library"
       liblvec.a.OBJS     = $(OBJS)
       liblvec.a.TYPE     = LIB
       liblvec.a.CP       = ../lib/liblvec.a

    #======================================================================
    #   Target: test_ivg                                                  |
    #----------------------------------------------------------------------

       test_ivg.NAME     = "Inheritance Vector test Program"
       test_ivg.TYPE     = C++
       test_ivg.OBJS     = test_ivg.o test_ivg_input.o
       test_ivg.DEP      = liblvec.a
       test_ivg.LDFLAGS  = -L../lib
       test_ivg.LDLIBS   = $(LIB_ALL)

    #======================================================================
    #   Target: test_nodes                                                |
    #----------------------------------------------------------------------

       test_nodes.NAME     = "Node Testing Program"
       test_nodes.TYPE     = C++
       test_nodes.OBJS     = test_nodes.o
       test_nodes.DEP      =
       test_nodes.LDFLAGS  = -L../lib


include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

