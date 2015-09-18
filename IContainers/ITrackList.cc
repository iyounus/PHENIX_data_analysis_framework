#include "ITrackList.hh"

#include "TClonesArray.h"
#include "ITrack.hh"


ClassImp(ITrackList)

#define NTRACKS 200

ITrackList::ITrackList()
{
  nTracks = 0;
  trackList = new TClonesArray("ITrack", NTRACKS);
}
//-----------------------------------------------------------------------------


ITrackList::~ITrackList()
{
  trackList->Clear();
  delete trackList;
}
//-----------------------------------------------------------------------------


void ITrackList::Reset()
{
  trackList->Clear();
  nTracks = 0;
}
//-----------------------------------------------------------------------------


ITrack* ITrackList::AddTrack()
{
  TClonesArray &tlist = *trackList;
  ITrack *track = new(tlist[nTracks++]) ITrack();
  return track;
}
//-----------------------------------------------------------------------------


ITrack* ITrackList::GetTrack(const unsigned int itrk)
{
  return (ITrack*)trackList->UncheckedAt(itrk);
}
//-----------------------------------------------------------------------------
