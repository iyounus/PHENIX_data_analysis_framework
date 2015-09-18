AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libIContainers.la

noinst_HEADERS = \
  IContainersLinkDef.h

libIContainers_la_SOURCES = \
  IHeader.cc       \
  IHeaderList.cc   \
  IPhoton.cc       \
  IPhotonList.cc   \
  ITrack.cc        \
  ITrackList.cc    \
  IMuon.cc         \
  IMuonList.cc     \
  IDiMuon.cc       \
  IDiMuonList.cc   \
  IElectron.cc     \
  IElectronList.cc \
  IBbc.cc          \
  IBbcList.cc      \
  IZdcSmd.cc       \
  IZdcSmdList.cc   \
  IContainers_Dict.C

libIContainers_la_LIBADD = \
  -L$(libdir)


# Rule for generating CINT dictionaries from class headers.
IContainers_Dict.C:   \
  IHeader.hh       \
  IHeaderList.hh   \
  IPhoton.hh       \
  IPhotonList.hh   \
  ITrack.hh        \
  ITrackList.hh    \
  IMuon.hh         \
  IMuonList.hh     \
  IDiMuon.hh       \
  IDiMuonList.hh   \
  IElectron.hh     \
  IElectronList.hh \
  IBbc.hh          \
  IBbcList.hh      \
  IZdcSmd.hh       \
  IZdcSmdList.hh   \
  IContainersLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
