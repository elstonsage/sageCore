# Makefile for Containers
#
#  Copyright (c) 2004  R.C. Elston

#==========================================================================
# General target options                                                  |
#--------------------------------------------------------------------------
#

  LSFROOT=..
  SAGEROOT=../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = Containers Library
  TARGET      =
  TARGETS     =
  TESTTARGETS = test_indexed_map test_sequence_map test_anyvector # test_process_mgr # test_functor
  VERSION     = 1.0
  TESTS       = true


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  DEP_SRCS = test_indexed_map.cpp test_anyvector.cpp test_sequence_map.cpp \
             test_process_mgr.cpp test_functor.cpp

    #======================================================================
    #   Target: test_indexed_map                                          |
    #----------------------------------------------------------------------

       test_indexed_map.NAME     = "Indexed Map Test"
       test_indexed_map.TYPE     = C++
       test_indexed_map.OBJS     = test_indexed_map.o
       test_indexed_map.CXXFLAGS = -L.

    #======================================================================
    #   Target: test_anyvector                                            |
    #----------------------------------------------------------------------

       test_anyvector.NAME     = "AnyVector test"
       test_anyvector.TYPE     = C++
       test_anyvector.OBJS     = test_anyvector.o

    #======================================================================
    #   Target: test_sequence_map                                         |
    #----------------------------------------------------------------------

       test_sequence_map.NAME     = "Indexed Map Test"
       test_sequence_map.TYPE     = C++
       test_sequence_map.OBJS     = test_sequence_map.o
       test_sequence_map.CXXFLAGS = -L.

    #======================================================================
    #   Target: test_process_mgr                                          |
    #----------------------------------------------------------------------

       test_process_mgr.NAME     = "ProcessMgr test"
       test_process_mgr.TYPE     = C++
       test_process_mgr.OBJS     = test_process_mgr.o
       test_process_mgr.CXXFLAGS = -L.

    #======================================================================
    #   Target: test_functor                                              |
    #----------------------------------------------------------------------

       test_functor.NAME     = "Functor / sequence test"
       test_functor.TYPE     = C++
       test_functor.OBJS     = test_functor.o
       test_functor.CXXFLAGS = -L.

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

