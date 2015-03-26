#==========================================================================
# Global Compiler options                                                 |
#--------------------------------------------------------------------------
#

  CC          = gcc
  CXX         = g++
  CFLAGS      = -DDEBUG_VERBOSE -I.
  CXXFLAGS    = -DDEBUG_VERBOSE -I. \
                -DNEEDS_LOG1P -DNEEDS_EXPM1 -DMINGW -DWINNT -DWINDOWS   \
                -DBROKEN_TEMPLATE_FRIENDS -DNEEDS_ERRORSTREAM_INIT      \
                -DBROKEN_VECTOR_ALLOCATOR -DNO_POSIX_SIGNALS            \
                -DUSE_PYTHON_GETOPT                                     \
                -ftemplate-depth-50                                     \
                -frtti -fexceptions                                     \
                -I../include                                            \
                -I$(SAGEROOT)/$(EXTERNDIR)/include              \
                -I$(SAGEROOT)/$(EXTERNDIR)/include/python1.5    \
                -I$(SAGEROOT)/$(EXTERNDIR)/include/boost_1_33_1

  LDLIBS      =
  LDFLAGS     = -DWINNT -DWINDOWS -L. -L$(SAGEROOT)/c++/lib -L$(SAGEROOT)/$(EXTERNDIR)/archive/$(ARCH)

  VERBOSE     =3
  OPTS        = sub
  MAKE        = make 
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

  EXE=.exe
  OBJ=.o
  LIB=.a
  LAPACK = -llapack -lblas -lg2c

  LAPACK_CFLAGS    = -g -I../include $(CFLAGS)

  XMLPP = -lxml++ -lxml2 -liconv -lz -L/usr/lib/w32api -lwsock32

  ERROR = -lerror -lerror_reg $(XMLPP)

  PYTHON = -lpython1.5 # -L$(SAGEROOT)/c++/experimental/python-embed/ -lpython1.5-mingw

  COVERAGE.CFLAGS   = -a -Wall
  COVERAGE.CXXFLAGS = -a -Wall
  COVERAGE.LDFLAGS  = -a

  DEBUG.CFLAGS      = -g -Wall
  DEBUG.CXXFLAGS    = -g -Wall

  STANDARD.CFLAGS   = -O2 -Wall
  STANDARD.CXXFLAGS = -O2 -Wall

  SPEED.CFLAGS      = -O -Wall
  SPEED.CXXFLAGS    = -O -Wall
  SPEED.LDFLAGS     = 

  RELEASE.CFLAGS    = -O2 -Wall -fomit-frame-pointer
  RELEASE.CXXFLAGS  = -O2 -Wall -fomit-frame-pointer
  RELEASE.LDFLAGS   = -static
  RELEASE.LDLIBS    = 
