#==========================================================================
# Global Compiler options                                                 |
#--------------------------------------------------------------------------
#

  CC            = /usr/bin/gcc
  CXX           = /usr/bin/g++
  CFLAGS        = -DDEBUG_VERBOSE -I.
  CXXFLAGS      = -DDEBUG_VERBOSE -I.
  LDLIBS        =  
  VERBOSE       =3
  OPTS          = sub
  MAKE          = make -j4
  MAKEDEPEND    = /usr/bin/X11/makedepend -Y. -I../lib -I/usr/include
  RM            = rm
  RM_RF         = -Rf
  CP            = cp
  MV            = mv
  AR            = ar
  AR_CREATE     = ruc
  AR_CXX_CREATE = ruc
  RANLIB        = true

  EXE=
  OBJ=.o
  LIB=.a
  LAPACK=-llapack -lblas -lf2c

  LAPACK_CFLAGS    = -g -I../include $(CFLAGS)

  XMLPP  = -lxml++ -lxml2 # -liconv -lpthread -lz -lm  

  PYTHON = -lpython1.5

  COVERAGE.CFLAGS      = -a -Wall
  COVERAGE.CXXFLAGS    = -a -Wall
  COVERAGE.LDFLAGS     = -a

  DEBUG.CFLAGS      = -g -Wall
  DEBUG.CXXFLAGS    = -g -Wall

  STANDARD.CFLAGS   = -O2 -Wall
  STANDARD.CXXFLAGS = -O2 -Wall

  SPEED.CFLAGS    = -O -Wall
  SPEED.CXXFLAGS  = -O -Wall

  RELEASE.CFLAGS    = -O2 -Wall 
  RELEASE.CXXFLAGS  = -O2 -Wall

#==========================================================================
# Global Compiler options for GCC on Mac/PowerPC                          |
#--------------------------------------------------------------------------

  CXXFLAGS   := $(CXXFLAGS)  -I/usr/include/                  \
               -I$(SAGEROOT)/c++/include -I../include  -Wall  \
               -DNEEDS_FORTRAN_MAIN -D_STL_NO_CONCEPT_CHECKS  \
               -I/usr/local/include                           \
               -I$(SAGEROOT)/$(EXTERNDIR)/include              \
               -I$(SAGEROOT)/$(EXTERNDIR)/include/python1.5    \
               -I$(SAGEROOT)/$(EXTERNDIR)/include/boost_1_49_0

  LDFLAGS    := -L. -L$(SAGEROOT)/c++/lib -L../lib -L$(SAGEROOT)/$(EXTERNDIR)/archive/$(ARCH) $(LDFLAGS)

  LINTFLAGS    := -i$(SAGEROOT)/c++/include -i../include          \
                  -dNEEDS_LOG1P -dNEEDS_EXPM1                     \
                  -i$(SAGEROOT)/$(EXTERNDIR)/include/python1.5     \
                  -i$(SAGEROOT)/$(EXTERNDIR)/include/boost_1_49_0  \
                  -d__i386__ -d__i386 -d__darwin__ -di386=1

  LINT_SRCS  = $($(TARGET).SRCS) $(SRCS)
  LINT_HEAD  = ${LINT_SRCS:.cpp=.h,} ${LINT_SRCS:.cpp=.ipp,}


  LINTINCLUDES :=  $(SAGEROOT)/config/lint/common.lnt $(SAGEROOT)/config/lint/co-g++.lnt
