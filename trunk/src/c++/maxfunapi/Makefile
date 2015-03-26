  SAGEROOT= ../..
  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = Maxfun
  TARGET      =
  TARGETS     = libmaxfunapi.a
  TESTTARGETS = libmaxfunapi.a maxtest maxtest2

  VERSION     = 1.0
  TESTS       = runall maxfunapi

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS        = APIMaxFunction.cpp           \
                CovarianceMatrix.cpp         \
                DebugCfg.cpp                 \
                GlobalFunctions.cpp          \
                Maximizer.cpp                \
                OutputFormatter.cpp          \
                Parameter.cpp                \
                ParameterInput.cpp           \
                ParameterMgr.cpp             \
                Results.cpp                  \
                RunCfg.cpp                   \
                SequenceCfg.cpp              \
                Submodel.cpp                 \
                TransformationSubmodel.cpp   \
                DivideByParamTransformer.cpp \
                AdditiveParamCalculator.cpp

  OBJS        = ${SRCS:.cpp=.o}

  DEP_SRCS    = maxtest.cpp maxtest2.cpp

  HEADERS     = ${SRCS:.cpp=.h}

    #======================================================================
    #   Target: libmaxfunapi.a                                               |
    #----------------------------------------------------------------------

       libmaxfunapi.a.NAME     = Maxfun API
       libmaxfunapi.a.TYPE     = LIB
       libmaxfunapi.a.OBJS     = ${OBJS}
       libmaxfunapi.a.CP       = ../lib/libmaxfunapi.a

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


    #======================================================================
    #   Target: maxtest2                                                     |
    #----------------------------------------------------------------------

       maxtest2.NAME      = Maxfun API Test
       maxtest2.INSTALL   = yes
       maxtest2.TYPE      = C++
       maxtest2.SRCS      = maxtest2.cpp
       maxtest2.OBJS      = maxtest2.o
       maxtest2.DEP       = libmaxfunapi.a
       maxtest2.LDLIBS    = $(LIB_ALL) 


include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

