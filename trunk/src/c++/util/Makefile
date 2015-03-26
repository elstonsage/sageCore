# Makefile for util 
#
#  Author: Geoff Wedig
#
#  History:   1.0   gcw  Initial creation              2004, Aug
#
#  Copyright (c) 2004  R.C. Elston
#

# NOTE:  At present (2004, Aug), there are no cpp files associated with the
#        S.A.G.E. util directory and it doesn't need building.  Despite
#        that, we need util to be part of the make system to allow for
#        testing, documentation and other utilities.  Thus we have this
#        directory and this empty Makefile. It is possible that in the
#        future, there will be need for this Makefile directly as other
#        utilities get added to the S.A.G.E. source.

#==========================================================================
# General target options                                                  |
#--------------------------------------------------------------------------
  LSFROOT=..
  SAGEROOT=../..

  include $(SAGEROOT)/config/Global.make

  TARGET_NAME = Miscellaneous S.A.G.E. Utilities
  TARGET      =
  TARGETS     = libutil.a
  TESTTARGETS = libutil.a test_disambiguator$(EXE) test_regex$(EXE) \
                test_autotrace$(EXE) \
                test_stringutils$(EXE) test_typeinfo$(EXE) \
                test_outline$(EXE)
  VERSION     = 
  TARPREFIX   = 
  TESTS       = runall util

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------  

  SRCS = AutoTrace.cpp

  DEP_SRCS = test_regx.cpp test_xmlparser.cpp test_stringutils.cpp \
             test_outline.cpp test_typeinfo.cpp test_objtracker.cpp \
             test_disambiguator.cpp

  #====================================================================== 
  #   Target: libutil.a                                                 |
  #---------------------------------------------------------------------- 

    libutil.a.NAME     = "Util library"
    libutil.a.TYPE     = LIB
    libutil.a.OBJS     = AutoTrace.o
    libutil.a.CP       = ../lib/libutil.a

  #====================================================================== 
  #   Target: test_regex                                                |
  #---------------------------------------------------------------------- 

    test_regex$(EXE).NAME          = Test regex functions
    test_regex$(EXE).INSTALL       = yes
    test_regex$(EXE).TYPE          = C++
    test_regex$(EXE).OBJS          = test_regex.o
    test_regex$(EXE).LDLIBS        = $(LIB_TOOLS)

  #====================================================================== 
  #   Target: test_xmlparser                                            |
  #---------------------------------------------------------------------- 

    test_xmlparser$(EXE).NAME          = Test xml parsing
    test_xmlparser$(EXE).INSTALL       = yes
    test_xmlparser$(EXE).TYPE          = C++
    test_xmlparser$(EXE).OBJS          = test_xmlparser.o
    test_xmlparser$(EXE).LDLIBS        = $(LIB_TOOLS)

  #====================================================================== 
  #   Target: test_stringutils                                          |
  #---------------------------------------------------------------------- 

    test_stringutils$(EXE).NAME          = Test string utils
    test_stringutils$(EXE).INSTALL       = yes
    test_stringutils$(EXE).TYPE          = C++
    test_stringutils$(EXE).OBJS          = test_stringutils.o
    test_stringutils$(EXE).LDLIBS        = $(LIB_TOOLS)

  #====================================================================== 
  #   Target: test_outline                                              |
  #---------------------------------------------------------------------- 

    test_outline$(EXE).NAME          = Test outline outline
    test_outline$(EXE).INSTALL       = yes
    test_outline$(EXE).TYPE          = C++
    test_outline$(EXE).OBJS          = test_outline.o
    test_outline$(EXE).LDLIBS        = $(LIB_TOOLS)

  #====================================================================== 
  #   Target: test_typeinfo                                             |
  #---------------------------------------------------------------------- 

    test_typeinfo$(EXE).NAME          = Test type info
    test_typeinfo$(EXE).INSTALL       = yes
    test_typeinfo$(EXE).TYPE          = C++
    test_typeinfo$(EXE).OBJS          = test_typeinfo.o
    test_typeinfo$(EXE).LDLIBS        = $(LIB_TOOLS)

  #====================================================================== 
  #   Target: test_objtracker                                           |
  #---------------------------------------------------------------------- 

    test_objtracker$(EXE).NAME          = Test objtracker
    test_objtracker$(EXE).INSTALL       = yes
    test_objtracker$(EXE).TYPE          = C++
    test_objtracker$(EXE).OBJS          = test_objtracker.o
    test_objtracker$(EXE).LDLIBS        = $(LIB_TOOLS)

  #====================================================================== 
  #   Target: test_disambiguator                                        |
  #---------------------------------------------------------------------- 

    test_disambiguator$(EXE).NAME          = Test disambiguator
    test_disambiguator$(EXE).INSTALL       = yes
    test_disambiguator$(EXE).TYPE          = C++
    test_disambiguator$(EXE).OBJS          = test_disambiguator.o
    test_disambiguator$(EXE).LDLIBS        =

  #======================================================================
  #   Target: test_autotrace                                            |
  #----------------------------------------------------------------------

    test_autotrace$(EXE).NAME          = Test auto tracer
    test_autotrace$(EXE).INSTALL       = yes
    test_autotrace$(EXE).CXXFLAGS      = -DAUTOTRACE_ON
    test_autotrace$(EXE).TYPE          = C++
    test_autotrace$(EXE).DEP           = libutil.a
    test_autotrace$(EXE).OBJS          = test_autotrace.o
    test_autotrace$(EXE).LDLIBS        = $(LIB_TOOLS)

include $(SAGEROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

