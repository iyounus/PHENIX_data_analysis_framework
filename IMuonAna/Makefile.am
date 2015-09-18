AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IContainers \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IPdstAna \
  -I$(includedir)  -I$(ROOTSYS)/include


lib_LTLIBRARIES = \
  libIMuonAna.la

noinst_HEADERS = \
  IMuonAnaLinkDef.h

libIMuonAna_la_SOURCES = \
  IMuonAna.cc            \
  IDiMuonAna.cc          \
  IDiMuonB2BAna.cc       \
  IMuonAna_Dict.C

libIMuonAna_la_LIBADD = \
  -L$(libdir)


# Rule for generating CINT dictionaries from class headers.
IMuonAna_Dict.C:    \
  IMuonAna.hh       \
  IDiMuonAna.hh     \
  IDiMuonB2BAna.hh  \
  IMuonAnaLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
