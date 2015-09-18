AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(HOME)/workdir/analysis/unm_pdst_ana/IContainers \
  -I$(includedir)  -I$(ROOTSYS)/include


lib_LTLIBRARIES = \
  libIPdstAna.la

noinst_HEADERS = \
  IPdstAnaLinkDef.h

libIPdstAna_la_SOURCES = \
  IPdst.cc        \
  IEvent.cc       \
  IConsts.cc      \
  IFillLookup.cc  \
  ISpinPattern.cc \
  IPdstAna_Dict.C

libIPdstAna_la_LIBADD = \
  -L$(libdir)


# Rule for generating CINT dictionaries from class headers.
IPdstAna_Dict.C:  \
  IPdst.hh        \
  IEvent.hh       \
  IConsts.hh      \
  IFillLookup.hh  \
  ISpinPattern.hh \
  IPdstAnaLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
