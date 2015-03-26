# Makefile for Sibpal
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

  TARGET_NAME = Sibpal
  TARGET      =
  TARGETS     = libsibpal.a sibpal$(EXE)
  TESTTARGETS = libsibpal.a sibpal$(EXE)
  VERSION     = 1.0
  TARPREFIX   = SIBPAL
  TESTS       = runall sibpal

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = definitions.h        sib_cluster.h        sib_matrix.h     \
                gls3.h               meantest_params.h    regress_params.h \
                parser.h             meantest.h           meantest_out.h   \
                regress_result.h     regress.h            regress_variants.h \
                regress_out.h        input.h              analysis.h       \
                sibpal.h

  SRCS        = sib_cluster.cpp      sib_matrix.cpp       gls3.cpp          \
                meantest_params.cpp  regress_params.cpp   parser.cpp        \
                meantest.cpp         meantest_out.cpp                       \
                regress_result.cpp   regress.cpp          regress_variants.cpp \
                regress_out.cpp      input.cpp            analysis.cpp

  DEP_SRCS    = sibpal.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libsibpal.a                                               |
    #----------------------------------------------------------------------

       libsibpal.a.NAME     = "Sibpal Library"
       libsibpal.a.TYPE     = LIB
       libsibpal.a.OBJS     = $(OBJS)
       libsibpal.a.CP       = ../lib/libsibpal.a

    #======================================================================
    #   Target: sibpal                                                    |
    #----------------------------------------------------------------------

       sibpal$(EXE).NAME      = Sibpal
       sibpal$(EXE).INSTALL   = yes
       sibpal$(EXE).TYPE      = C++
       sibpal$(EXE).OBJS      = sibpal.o
       sibpal$(EXE).DEP       = libsibpal.a
       sibpal$(EXE).LDLIBS    = -lsibpal -lpalbase $(LIB_ALL)

include $(SAGEROOT)/config/Rules.make

test_update:
	-@for i in test?; do                   \
            cd $$i                           ; \
            for j in *.out; do                 \
              mv $$j `basename $$j .out`.exp ; \
            done                             ; \
            cd ..                            ; \
          done

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

