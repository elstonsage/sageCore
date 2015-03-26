#==========================================================================
# General target options                                                  |
#--------------------------------------------------------------------------
#

  LSFROOT=..
  SAGEROOT=../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = Functors Library
  TARGET      =
  TARGETS     =
  TESTTARGETS = test_functor
  VERSION     = 1.0
  TESTS       = runall functors


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  DEP_SRCS    = test_functor.cpp

  OBJS        = ${SRCS:.cpp=.o}

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

