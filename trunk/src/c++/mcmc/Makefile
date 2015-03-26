# Makefile for MCMC basic objects

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

  TARGET_NAME = "MCMC basic objects"
  TARGET      =
  TARGETS     = libmcmc.a
  TESTTARGETS = libmcmc.a test_mcmc$(EXE)
  VERSION     = 1.1
  TESTS       = runall mcmc

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     =  mcmc_meiosis_map.h     \
                 mcmc_params.h          \
                 mcmc_parser.h          \
                 common_allele_set.h    \
                 mcmc_data_accessor.h   \
                 founder_allele_graph.h \
                 marker_likelihoods.h   \
                 recomb_calculator.h    \
                 starting_state.h       \
                 mcmc_simulator.h      

  SRCS        = ${HEADERS:.h=.cpp }

  DEP_SRCS    = test_common_allele_set.cpp founder_allele_graph.cpp \
                test_mcmc.cpp test_mcmc_params.cpp test_mcmc_input.cpp \
                test_mcmc_parser.cpp test_mcmc_analysis.cpp

  OBJS        = ${SRCS:.cpp=.o}


    #======================================================================
    #   Target: libmcmc.a                                                 |
    #----------------------------------------------------------------------

       libmcmc.a.NAME     = "MCMC simulation Library"
       libmcmc.a.TYPE     = LIB
       libmcmc.a.OBJS     = $(OBJS)
       libmcmc.a.CP       = ../lib/libmcmc.a

    #======================================================================
    #   Target: test_common_allele_set                                                 |
    #----------------------------------------------------------------------

       test_common_allele_set$(EXE).NAME     = "Test the CommonAlleleSet component"
       test_common_allele_set$(EXE).INSTALL  = yes
       test_common_allele_set$(EXE).TYPE     = C++
       test_common_allele_set$(EXE).OBJS     = test_common_allele_set.o founder_allele_graph.o
       test_common_allele_set$(EXE).LDLIBS   = $(LIB_ALL)

    #======================================================================
    #   Target: test_mcmc                                                 |
    #----------------------------------------------------------------------

       test_mcmc$(EXE).NAME     = "MCMC Components Testing Program"
       test_mcmc$(EXE).INSTALL  = yes
       test_mcmc$(EXE).TYPE     = C++
       test_mcmc$(EXE).DEP      = libmcmc.a
       test_mcmc$(EXE).OBJS     = test_mcmc.o test_mcmc_params.o test_mcmc_input.o \
                                  test_mcmc_parser.o test_mcmc_analysis.o
       test_mcmc$(EXE).LDLIBS   = $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

