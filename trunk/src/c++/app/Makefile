# Makefile for Pair Generator Main Body
#
#  Authors: Kevin Jacobs (jacobs@darwin.cwru.edu)
#           Geoff Wedig  (wedig@darwin.cwru.edu)
#
#  History:   1.0   kbj  Initial implementation        Apr  5 1996
#             1.01  kbj  Refinements and confusion     Jun 10 1996
#             1.1   kbj  Finally got it right          Jul  1 1996
#             1.11  gcw  Relative Pair port            Jul 17 1996
#             1.2   gcw  Added Global.make, Rules.make Aug 19 1996
#             1.21  gcw  Written for RPped             Aug 20 1996
#             1.3   gcw  New Source handling           Jan 06 1997
#
#  Copyright (c) 1996  R.C. Elston
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

  SAGEROOT= ../..
  include $(SAGEROOT)/config/Global.make

  TTARGETS    = 

  TARGET_NAME = Application Framework Code
  TARGET      =
  TARGETS     = libapp.a 
  TESTTARGETS = libapp.a test_app
  VERSION     = 1.0
  TESTS       = runall app

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS = SAGEapp.cpp SAGEapp_version_bank.cpp output_streams.cpp aparser.cpp VersionNumber.cpp
  OBJS = ${SRCS:%.cpp=%.o}

  DEP_SRCS = test_app.cpp

  #========================================================================
  # APP Library Sources                                                   |
  #------------------------------------------------------------------------

    #======================================================================
    #   Target: libapp.a                                                  |
    #----------------------------------------------------------------------

       libapp.a.NAME     = "Application Framework Code"
       libapp.a.OBJS     = $(OBJS)
       libapp.a.TYPE     = LIB
       libapp.a.CP       = ../lib/libapp.a


TEMPUNAME := $(shell uname -s)

EXTRALIB := -lboost_regex

ifeq ($(ARCH), i386-mingw-gcc)

  EXTRALIB := -L"C:/Dev-Cpp/lib" -liphlpapi

endif

ifeq ($(TEMPUNAME), SunOS)

  EXTRALIB := -lsocket -lnsl

endif

    #======================================================================
    #   Target: test_app
    #----------------------------------------------------------------------

       test_app.NAME          = Age of Onset Application
       test_app.INSTALL       = yes
       test_app.TYPE          = C++
       test_app.OBJS          = test_app.o
       test_app.DEP           = libapp.a
       test_app.LDLIBS        = $(LIB_DATA_CLEANING) $(EXTRALIB)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

