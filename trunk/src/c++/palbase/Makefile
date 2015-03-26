# Makefile for Bases for Pair Analyses (SIBPAL, LODPAL, RELPAL)
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

  TARGET_NAME = "Bases for Pair Analyses"
  TARGET      =
  TARGETS     = libpalbase.a 
  TESTTARGETS = libpalbase.a #testrelpair$(EXE)
  TARPREFIX   = PALBASE
  TESTS       = true


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = util.h         pair_info.h       relative_pairs.h  \
                rel_pair.h     pal_ibd.h         pair_info_file.h  \
                pair_filter.h

  SRCS        = relative_pairs.cpp  pair_filter.cpp   pair_info_file.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libpalbase.a                                            |
    #----------------------------------------------------------------------

       libpalbase.a.NAME     = New Pairs Library
       libpalbase.a.TYPE     = LIB
       libpalbase.a.OBJS     = $(OBJS)
       libpalbase.a.CP       = ../lib/libpalbase.a


include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)

