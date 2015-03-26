#==========================================================================
# Setups for multi-builds.                                                |
#--------------------------------------------------------------------------

#BUILD_OBJS   := $(OBJS:%=$(BUILDDIR)/%)
TOJ           := $($(TARGET).OBJS)
TARGET_OBJS   := $(TOJ:%=$(BUILDDIR)/%)
GPATH          = $(BUILDDIR):.

# Include directive for platform-specific detection
include ${SAGEROOT}/config/PlatformDetection.make

#==========================================================================
# Primary rule.  Do not touch!!!!!!                                       |
#--------------------------------------------------------------------------

all: STATUS recurse

#==========================================================================
# Rule to generate blob source file.                                      |
#--------------------------------------------------------------------------

recurse:
	@if [ "$(TARGETS)" ]; then                                                             \
	  $(MAKE) --no-print-directory recurse2 AUTOTRACE=$(AUTOTRACE) BUILD=$(BUILD) $(OPTS); \
	fi

recurse2:
	@unset MAKEFLAGS;                                                                        \
	set -e;                                                                                  \
	for i in $(TARGETS) ; do                                                                 \
          $(MAKE) --no-print-directory TARGET=$$i AUTOTRACE=$(AUTOTRACE) BUILD=$(BUILD) $(OPTS); \
        done;                                                                                    \
	if [ $(VERBOSE) -gt 2 ]; then                                                            \
	  ls $(TARGETS) ;                                                                        \
	fi                                                                                       \

# Leaf rule
sub: LSTATUS $(TARGET)

# Main target rule
$(TARGET): $($(TARGET).DEP) $($(TARGET).LIBS) $(TARGET_OBJS)
	@BUILD_LDFLAGS="$(LDFLAGS) $($(BUILD).LDFLAGS) $($(TARGET).LDFLAGS)";      \
	BUILD_CXXFLAGS="$(CXXFLAGS) $($(BUILD).CXXFLAGS) $($(TARGET).CXXFLAGS)";   \
	BUILD_CFLAGS="$(CFLAGS) $($(BUILD).CFLAGS) $($(TARGET).CFLAGS)"        ; \
	BUILD_LDLIBS="$(LDLIBS) $($(BUILD).LDLIBS) $($(TARGET).LDLIBS)"        ; \
	if [ "$($(TARGET).TYPE)" != "LIB" -a "$($(TARGET).TYPE)" != "CLIB" ]; then  \
          if [ "$($(TARGET).TYPE)" = "C++" ]; then                          \
	    echo "Linking                 (C++) $@";                        \
	    $(CXX) $$BUILD_LDFLAGS -o $@ $(TARGET_OBJS) $$BUILD_LDLIBS ;    \
	  else                                                              \
	    echo "Linking                 (C)   $@" ;                       \
	    $(CC)  $$BUILD_LDFLAGS -o $@ $(TARGET_OBJS) $$BUILD_LDLIBS ;    \
	  fi ;                                                              \
	else                                                                \
           if [ "$($(TARGET).CLOSURE)" ]; then                              \
	     echo "Currently not supported.";                               \
	     exit 1;                                                        \
           fi ;                                                             \
	   printf "Making Library..." ;                                     \
	   LOCAL_DIR=`pwd`;                                                 \
	   cd $(BUILDDIR) ;                                                 \
	   if [ "$($(TARGET).TYPE)" = "CLIB" ]; then                        \
             $(AR) $(AR_CREATE) $(TARGET) $($(TARGET).OBJS) ;               \
             $(RANLIB) $(TARGET) ;                                          \
	   elif [ -n "$(AR_CXX)" ]; then                                    \
	     $(AR_CXX) $(AR_CXX_CREATE) $(TARGET) $($(TARGET).OBJS)         \
                                           $$BUILD_LDFLAGS $$BUILD_LDLIBS;  \
           else                                                             \
             $(AR) $(AR_CREATE) $(TARGET) $($(TARGET).OBJS) ;               \
             $(RANLIB) $(TARGET) ;                                          \
	   fi ;                                                             \
	   mv $(TARGET) $$LOCAL_DIR ;                                       \
	   printf " done.\n" ;                                              \
	   cd $$LOCAL_DIR ;                                                 \
	fi ;                                                                \
	if [ "$($(TARGET).TNAME)" ]; then                                   \
	  ln -fs $(TARGET) $($(TARGET).TNAME);                              \
	fi ;                                                                \
	if [ "$($(TARGET).CP)" ]; then                                      \
	  $(CP) $(TARGET) $($(TARGET).CP);                                  \
	fi ;                                                                \
	if [ "$($(TARGET).MV)" ]; then                                      \
	  $(MV) $(TARGET) $($(TARGET).MV);                                  \
	fi ;                                                                \
	if [ "$($(TARGET).RM)" ]; then                                      \
	  $(RM) $($(TARGET).RM) 2> /dev/null;                               \
	fi

install:
	@$(MAKE) OPTS=install-target

install-target: sub
	@if [ -n "$(INSTALL_DIR)" -a                          \
              -n "$($(TARGET).INSTALL)" ] ; then              \
          if [ '!' -d $(INSTALL_DIR) ]; then                  \
            mkdir $(INSTALL_DIR);                             \
          fi ;                                                \
	  echo "Installing $(TARGET) to $(INSTALL_DIR)...";   \
	  if [ \! -d $(INSTALL_DIR) ]; then                   \
	    mkdir $(INSTALL_DIR);                             \
	  fi ;                                                \
	  if [ $($(TARGET).TNAME) ] ; then                    \
	    cp $(TARGET) $(INSTALL_DIR)/$($(TARGET).TNAME);   \
	  else                                                \
	    cp $(TARGET) $(INSTALL_DIR)/$(TARGET);            \
	  fi;                                                 \
	fi

# Utility rules

strip:
	-@strip $(TARGETS)
	@ls -s $(TARGETS)

relink:
	@-rm $(TARGETS)
	@-$(MAKE) all

lint:
	@if [ "$($(TARGET).SRCS) $(SRCS)" ]; then                       \
	    echo \-libh\($(LINT_HEAD)\) >tmpfile~  ; \
	    sed -e 's/ / $(CURRENT_DIR)\//g' -e 's/.h, /.h) \-libh(/g'  \
	        -e 's/.ipp, /.ipp) -libh(/g' \
	        -e 's/,)/)/g' tmpfile~ >tmp2~.lnt ; \
	    $(FLINT) $(LINTINCLUDES) tmp2~.lnt $(LINTFLAGS) \
		$(LINT_SRCS)  ;                                         \
	    rm tmpfile~ tmp2~.lnt ; \
	fi

docs: docs2
	@echo Done
docs2:
	@../../../util/doc/moddocs.py

clean: $(CLEAN_TASKS) clean_target clean_source

cleanall: $(CLEANALL_TASKS) clean clean_target_repository

clean_source: $(CLEAN_TASKS)
	-@$(RM) $(TARGETS) $(TEST_TARGETS) $(TESTTARGETS) *.a *.bak *.exe *~ 2> /dev/null

clean_target:
	-@$(RM) $(BUILDDIR)/*$(OBJ) $(BUILDDIR)/*$(LIB)  $(BUILDDIR)/*.d 2> /dev/null
	-@$(RM) $(BUILDDIR)/KCC_files/*$(OBJ) $(BUILDDIR)/ti_files/*$(OBJ) 2> /dev/null

clean_target_repository:
	-@$(RM) $(BUILDDIR)/*.ii $(BUILDDIR)/KCC_files/*.ii 2>/dev/null
	-@$(RM) $(BUILDDIR)/ti_files/*.ti $(BUILDDIR)/ti_files/*.ii 2>/dev/null
	-@$(RM) $(BUILDDIR)/ti_files/*.ci 2>/dev/null

# New nicer rules to replace the builtins for pattern based dependencies

.SUFFIXES:
.SUFFIXES: .cc .cpp .c $(OBJ) .h

$(BUILDDIR)/%$(OBJ): %.cpp
	@echo "Compiling               (C++) $*.cpp"
	@BUILD_CXXFLAGS="$(CXXFLAGS) $($(BUILD).CXXFLAGS) $($(TARGET).CXXFLAGS)";                 \
	BUILD_LDFLAGS="$(LDFLAGS) $($(BUILD).LDFLAGS) $($(TARGET).LDFLAGS)"     ;                 \
	BUILD_CXXFLAGS="$(CXXFLAGS) $($(BUILD).CXXFLAGS) $($(TARGET).CXXFLAGS)" ;                 \
	BUILD_LDLIBS="$(LDLIBS) $($(BUILD).LDLIBS) $($(TARGET).LDLIBS)"         ;                 \
	$(CXX) -MD $$BUILD_CXXFLAGS  -o $(BUILDDIR)/$*$(OBJ) -c $*.cpp                                \
          -DBUILD=\""$$BUILD"\" -DCXXFLAGS=\""$$BUILD_CXXFLAGS"\" -DLDFLAGS=\""$$BUILD_LDFLAGS"\" \
          -DLDLIBS=\""$$BUILD_LDLIBS"\" -DBUILD_DATE=\""`date '+%d %b %Y'`"\"
	@sed -e 's/^ *//' -e '/\/usr\/.*/ d' < $(BUILDDIR)/$*.d > $(BUILDDIR)/$*.dtmp; \
	cp $(BUILDDIR)/$*.dtmp $(BUILDDIR)/$*.d; \
	echo "" >> $(BUILDDIR)/$*.d; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
            -e 's/^ *//' -e '/\/usr\/.*/ d' \
            -e '/^$$/ d' -e 's/$$/ :/' < $(BUILDDIR)/$*.dtmp >> $(BUILDDIR)/$*.d; \
	rm $(BUILDDIR)/$*.dtmp

$(BUILDDIR)/%$(OBJ): %.c
	@echo "Compiling               (C)   $*.c"
	@BUILD_CFLAGS="$(CFLAGS) $($(BUILD).CFLAGS) $($(TARGET).CFLAGS)"; \
	$(CC) $$BUILD_CFLAGS -c $*.c -o $(BUILDDIR)/$*$(OBJ)

# Print out nice status information during build

# Global status
STATUS:
	@if [ "$(VERBOSE)" -gt 2 ]; then                                    \
	  BUILD_LDFLAGS="$(LDFLAGS) $($(BUILD).LDFLAGS)";                   \
	  BUILD_CXXFLAGS="$(CXXFLAGS) $($(BUILD).CXXFLAGS)";                \
	  BUILD_CFLAGS="$(CFLAGS) $($(BUILD).CFLAGS)";                      \
	  BUILD_LDLIBS="$(LDLIBS) $($(BUILD).LDLIBS)";                      \
	  echo "";                                                          \
	  echo "Making $(TARGET_NAME) for $(ARCH) $(BUILD) in $(BUILDDIR)"; \
	  echo "CC       = $(CC)";                                          \
	  echo "CXX      = $(CXX)";                                         \
	  echo "CFLAGS   = $$BUILD_CFLAGS";                                 \
	  echo "CXXFLAGS = $$BUILD_CXXFLAGS";                               \
	  echo "LDFLAGS  = $$BUILD_LDFLAGS";                                \
	  echo "LDLIBS   = $$BUILD_LDLIBS";                                 \
	  echo "";                                                          \
	fi

# Local status

LSTATUS:
	@if [ '$(VERBOSE)' -gt 1 ]; then \
	  if [ '$(TARGET)' ]; then \
	    echo 'Making $(TARGET): $($(TARGET).NAME)'; fi;   \
	  if [ '$($(TARGET).CFLAGS)' ]; then \
	    echo '  CFLAGS   = $($(TARGET).CFLAGS)'   ; fi; \
	  if [ '$($(TARGET).CXXFLAGS)' ]; then \
	    echo '  CXXFLAGS = $($(TARGET).CXXFLAGS)' ; fi; \
	  if [ '$($(TARGET).LDFLAGS)' ]; then \
	    echo '  LDFLAGS  = $($(TARGET).LDFLAGS)'  ; fi; \
	  if [ '$($(TARGET).LDLIBS)' ]; then \
	    echo '  LDLIBS   = $($(TARGET).LDLIBS)'   ; fi; \
	  echo "" ; \
	fi
        
DEPEND:
	echo "DEPEND" $(TARGET)

# Testing rules

testing: test

test: testbuild $(TEST_TASKS)
	-@printf 'Testing $(TARGET_NAME): '
	-@set -e;                                  \
        if $(TESTS); then                          \
	  echo "Passed";                           \
	else                                       \
	  echo "FAILED!";                          \
	fi

testbuild:
	-@if [ -n "$(TESTTARGETS)" ]; then \
	  $(MAKE) TARGETS="$(TESTTARGETS)" BUILD=$(BUILD); \
	fi

help:
	@cat $(SAGEROOT)/config/Makefile.help
