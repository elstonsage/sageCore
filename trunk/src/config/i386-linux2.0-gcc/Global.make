#==========================================================================
# Global Compiler options                                                 |
#--------------------------------------------------------------------------
#
#

  CC          = gcc
  CXX         = g++
  CFLAGS      = -DDEBUG_VERBOSE -I.
  CXXFLAGS    = -DDEBUG_VERBOSE -I.
  LDLIBS      = # -lefence
  VERBOSE     =3
  OPTS        = sub
  MAKE        = make -j4
  MAKEDEPEND  = /usr/bin/X11/makedepend -Y. -I../lib -I/usr/include
  RM          = rm
  RM_RF       = -Rf
  CP          = cp
  MV          = mv
  AR          = ar
# AR_CXX      = ar
  AR_CREATE   = ruc
  AR_CXX_CREATE = ruc
  RANLIB      = true
  FLINT       = ~/FlexLint/flint

  EXE=
  OBJ=.o
  LIB=.a
  LAPACK=-llapack -lblas -lf2c

  EXCEPTIONS = --no_exceptions --no_rtti

  PYTHON = -lpython1.5 -ldl

  XMLPP  = -lxml++ -lxml2 -liconv -lpthread -lz -lm  

  FLEXLM =  -L/usr/local/lib/ /usr/local/lib/lm_new.o -llmgr -lcrvs -lsb

  LIB_PLATSPEC = -lpthread -lc

  COVERAGE.CFLAGS      = -a
  COVERAGE.CXXFLAGS    = -a
  COVERAGE.LDFLAGS     = -a

  PROFILE.CXXFLAGS = -O3 -pg -g
  PROFILE.LDFLAGS = -pg

  DEBUG.CFLAGS      = -g
  DEBUG.CXXFLAGS    = -g

  STANDARD.CFLAGS   = -O
  STANDARD.CXXFLAGS = -O

  RELEASE.CFLAGS    = -O3 -fomit-frame-pointer -funroll-loops
  RELEASE.CXXFLAGS  = -O3 -fomit-frame-pointer -funroll-loops
  RELEASE.LDFLAGS   = -static

#==========================================================================
# Global Compiler options for GCC on Linux/i386                           |
#--------------------------------------------------------------------------

  CXXFLAGS   := $(CXXFLAGS)  -I/usr/include/  \
               -I$(SAGEROOT)/c++/include -I../include  -Wall \
               -DNEEDS_FORTRAN_MAIN -D_STL_NO_CONCEPT_CHECKS      \
               -I/usr/local/include                           \
               -I$(SAGEROOT)/$(EXTERNDIR)/include              \
               -I$(SAGEROOT)/$(EXTERNDIR)/include/python1.5    \
               -I$(SAGEROOT)/$(EXTERNDIR)/include/boost_1_34_1

  LDFLAGS    :=-L. -L$(SAGEROOT)/c++/lib -L../lib -L$(SAGEROOT)/$(EXTERNDIR)/archive/$(ARCH) $(LDFLAGS)
  LDLIBS     := $(LDLIBS) -lieee

  LINTFLAGS    := -i$(SAGEROOT)/c++/include -i../include \
                  -dNEEDS_LOG1P -dNEEDS_EXPM1                      \
                  -i$(SAGEROOT)/$(EXTERNDIR)/include/python1.5    \
                  -i$(SAGEROOT)/$(EXTERNDIR)/include/boost_1_34_1   \
                  -d__i386__ -d__i386 -d__linux__ -di386=1

  LINT_SRCS  = $($(TARGET).SRCS) $(SRCS)
  LINT_HEAD = ${LINT_SRCS:.cpp=.h,} ${LINT_SRCS:.cpp=.ipp,}


  LINTINCLUDES :=  $(SAGEROOT)/config/lint/common.lnt $(SAGEROOT)/config/lint/co-g++.lnt
