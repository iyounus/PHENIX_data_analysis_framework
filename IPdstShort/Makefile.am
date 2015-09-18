AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IContainers \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IPdstAna \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IPdstShort \
  -I$(includedir) -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libIPdstShort.la

noinst_HEADERS = \
  IPdstShortLinkDef.h

libIPdstShort_la_SOURCES = \
  IPdstShort.cc            \
  IPdstShort_Dict.C

libIPdstShort_la_LIBADD = \
  -L$(libdir)

# Rule for generating CINT dictionaries from class headers.
IPdstShort_Dict.C: \
  IPdstShort.hh    \
  IPdstShortLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
