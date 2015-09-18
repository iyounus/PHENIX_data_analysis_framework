AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IContainers \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IPdstAna \
  -I$(includedir)  -I$(ROOTSYS)/include


lib_LTLIBRARIES = \
  libIGammaHadronAna.la

noinst_HEADERS = \
  IGammaHadronAnaLinkDef.h

libIGammaHadronAna_la_SOURCES = \
  IGammaHadron.cc           \
  IGammaHadronAna_Dict.C

libIGammaHadronAna_la_LIBADD = \
  -L$(libdir)


# Rule for generating CINT dictionaries from class headers.
IGammaHadronAna_Dict.C:   \
  IGammaHadron.hh     \
  IGammaHadronAnaLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
