AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IContainers \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IPdstAna \
  -I$(includedir)  -I$(ROOTSYS)/include


lib_LTLIBRARIES = \
  libIPi0Ana.la

noinst_HEADERS = \
  IPi0AnaLinkDef.h

libIPi0Ana_la_SOURCES = \
  IPi0Reco.cc           \
  IPi0Ana_Dict.C

libIPi0Ana_la_LIBADD = \
  -L$(libdir)


# Rule for generating CINT dictionaries from class headers.
IPi0Ana_Dict.C:   \
  IPi0Reco.hh     \
  IPi0AnaLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
