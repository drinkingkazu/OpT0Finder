
ifndef FMATCH_BASEDIR
ERROR_MESSAGE := $(error FMATCH_BASEDIR is not set... run configure.sh!)
endif

OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(FMATCH_BASEDIR)/Makefile/Makefile.${OSNAME}

SUBDIRS := Base PhotonLibrary #DataFormat Processor CPPUtil
#ifeq ($(LARCV_NUMPY),1)
#CORE_SUBDIRS += PyUtil
#endif

.phony: all clean

all: obj lib
#all: obj

clean: clean_app clean_core
	@rm -f $(FMATCH_LIBDIR)/libflashmatch.so
clean_core:
	@for i in $(SUBDIRS); do ( echo "" && echo "Cleaning $$i..." && cd $(FMATCH_DIR)/$$i && rm -rf $(FMATCH_BUILDDIR)/$$i && rm -rf $(FMATCH_BUILDDIR)/lib/*$ii.* ) || exit $$?; done

obj:
	@echo
	@echo Building ...
	@echo
	@for i in $(SUBDIRS); do ( echo "" && echo "Compiling $$i..." && cd $(FMATCH_DIR)/$$i && $(MAKE) ) || exit $$?; done

lib: obj
	@echo
	@if [ `${FMATCH_PYTHON} ${FMATCH_BASEDIR}/bin/libarg.py build` ]; then \
	   echo Linking library...; \
	   $(SOMAKER) $(SOFLAGS) $(shell ${FMATCH_PYTHON} $(FMATCH_BASEDIR)/bin/libarg.py); \
	 else \
	   echo Nothing to be done for lib...; \
	fi
