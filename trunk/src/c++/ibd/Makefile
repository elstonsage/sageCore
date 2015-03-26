# Makefile for IBD Objects
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

  TARGET_NAME = "Identical By Descent"
  TARGETS     = libibd.a
  TESTTARGETS = libibd.a
  VERSION     = 1.0
  TARPREFIX   = IBD
  TESTS       = true

#==========================================================================
# Source/object file lists                                                |

#--------------------------------------------------------------------------

  HEADERS     = ibd.h prior_ibd.h   ibdfile.h            ibd_analysis.h \
                basic_storage_ibd.h exact_ibd_analysis.h

  SRCS        = prior_ibd.cpp       ibdfile.cpp          ibd_analysis.cpp \
                exact_ibd_analysis.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libibd.a                                                  |
    #----------------------------------------------------------------------

       libibd.a.NAME     = "Identical By Descent Library"
       libibd.a.OBJS     = $(OBJS)
       libibd.a.TYPE     = LIB
       libibd.a.CP       = ../lib/libibd.a

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)

