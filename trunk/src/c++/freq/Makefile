#==========================================================================
# General target options                                                  |
#--------------------------------------------------------------------------

  SAGEROOT=../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME         = FREQ
  TARGET              =
  TARGETS             = freq$(EXE)
  TESTTARGETS         = freq$(EXE)
  TESTS               = runall freq

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS                = Sample.cpp Parser.cpp Output.cpp Peeler.cpp \
                        SingletonCalc.cpp LikelihoodCalculator.cpp \
                        Estimator.cpp freq.cpp

  DEP_SRCS            = freq.cpp

  OBJS                = ${SRCS:.cpp=.o}

  #======================================================================
  #   Target: libdecipher.a                                             |
  #----------------------------------------------------------------------

    freq$(EXE).NAME     = "FREQ"
    freq$(EXE).INSTALL  = yes
    freq$(EXE).TYPE     = C++
    freq$(EXE).OBJS     = $(OBJS)
    freq$(EXE).LDLIBS   = $(LIB_ALL)  # -lefence

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

