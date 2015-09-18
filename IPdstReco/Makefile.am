AUTOMAKE_OPTIONS = foreign

# I'm trying to include IContainers in two different ways. I'm not sure which one of these works for build env.
INCLUDES = \
  -I$(HOME)/workdir/analysis/unm_pdst/IContainers \
  -I$(includedir) -I$(srcdir)/../IContainers -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libIPdstReco.la 

noinst_HEADERS = \
  IPdstRecoLinkDef.h


libIPdstReco_la_SOURCES = \
  IPdstReco.cc \
  IPdstReco_Dict.C

libIPdstReco_la_LIBADD = \
  -L$(libdir)
  #-lfun4allfuncs


# Rule for generating CINT dictionaries from class headers.
IPdstReco_Dict.C: \
  IPdstReco.hh \
  IPdstRecoLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
