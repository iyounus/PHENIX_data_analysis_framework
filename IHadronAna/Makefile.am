AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IContainers \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IPdstAna \
  -I$(includedir)  -I$(ROOTSYS)/include


lib_LTLIBRARIES = \
  libIHadronAna.la

noinst_HEADERS = \
  IHadronAnaLinkDef.h

libIHadronAna_la_SOURCES = \
  IHadronAna.cc           \
  IHadronAna_Dict.C

libIHadronAna_la_LIBADD = \
  -L$(libdir)


# Rule for generating CINT dictionaries from class headers.
IHadronAna_Dict.C:   \
  IHadronAna.hh     \
  IHadronAnaLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
