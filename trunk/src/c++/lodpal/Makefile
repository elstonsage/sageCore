# Makefile for LODPAL
#
#  Authors: Kevin Jacobs (jacobs@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Oct 22 1998
#
#  Copyright (c) 1998  R.C. Elston
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
  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = LODSCORE test of Affected Relative Pair
  TARGET      =
  TARGETS     = liblodpal.a lodpal$(EXE)
  TESTTARGETS = liblodpal.a lodpal$(EXE) #test_maxfunapi
  TARPREFIX   = LODPAL
  TESTS       = runall lodpal


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = lodpal_parser.h         lodpal_params.h         lodpal_pairs.h \
                lodpal_input.h          arp_base_analysis.h     arp_one_analysis.h \
                arp_two_analysis.h      arp_po_one_analysis.h   arp_po_two_analysis.h \
                arp_x_one_analysis.h    arp_x_two_analysis.h    dsp_one_analysis.h \
                lodpal_result.h         lodpal_out.h            lodpal_analysis.h \
                lodpal.h

  SRCS        = lodpal_parser.cpp       lodpal_params.cpp       lodpal_pairs.cpp \
                lodpal_input.cpp        arp_base_analysis.cpp   arp_one_analysis.cpp \
                arp_two_analysis.cpp    arp_po_one_analysis.cpp arp_po_two_analysis.cpp \
                arp_x_one_analysis.cpp  arp_x_two_analysis.cpp  dsp_one_analysis.cpp \
                lodpal_out_text.cpp     lodpal_out_diag.cpp     lodpal_out_csv.cpp \
                lodpal_out_x.cpp        lodpal_analysis.cpp

  DEP_SRCS    = lodpal.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: liblodpal.a                                               |
    #----------------------------------------------------------------------

       liblodpal.a.NAME     = "Lodpal Library"
       liblodpal.a.TYPE     = LIB
       liblodpal.a.OBJS     = $(OBJS)
       liblodpal.a.CP       = ../lib/liblodpal.a

    #======================================================================
    #   Target: lodpal                                                    |
    #----------------------------------------------------------------------

       lodpal$(EXE).NAME      = "LODSCORE test of Affected Relative Pair"
       lodpal$(EXE).TYPE      = C++
       lodpal$(EXE).INSTALL   = yes
       lodpal$(EXE).CXXFLAGS  =
       lodpal$(EXE).OBJS      = lodpal.o
       lodpal$(EXE).LDLIBS    = -llodpal -lpalbase $(LIB_ALL) 

    #======================================================================
    #   Target: test_maxfunapi                                            |
    #----------------------------------------------------------------------

       test_maxfunapi.NAME      = "Maxfun API Test"
       test_maxfunapi.TYPE      = C++
       test_maxfunapi.INSTALL   = yes
       test_maxfunapi.CXXFLAGS  =
       test_maxfunapi.OBJS      = test_maxfunapi.o
       test_maxfunapi.LDLIBS    = $(LIB_ALL) 


include $(SAGEROOT)/config/Rules.make



# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

