//-----------------------------------------------------------
// IHeader.hh
//       Created April 11 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef ITRACKLIST_HH
#define ITRACKLIST_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TClonesArray;
class ITrack;

class ITrackList : public TObject
{
public:
  ITrackList();
  virtual ~ITrackList();

  void Reset();

  unsigned short GetNTracks() const { return nTracks; }

  ITrack* AddTrack();
  ITrack* GetTrack(const unsigned int itrack);

protected:
  unsigned short nTracks;
  TClonesArray *trackList;

private:
  ClassDef(ITrackList,1);
};
#endif
