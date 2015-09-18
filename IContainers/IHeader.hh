//-----------------------------------------------------------
// IHeader.hh
//       Created April 11 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IHEADER_HH
#define IHEADER_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class IHeader : public TObject
{
private:
  int   RunID;
  int   EventID;
  float VtxZ;

  // spin information
  short GL1CrossingID;
  short SpinGL1CrossingID;  // corrected crossing ID for bunch shift

  unsigned int TrigScaled;
  unsigned int TrigLive;


public:
  IHeader(){;}
  ~IHeader(){;}


  int   GetRunID()   const {return RunID;}  
  int   GetEventID() const {return EventID;} 
  float GetZVertex() const {return VtxZ;}

  short GetGL1CrossingID()     const {return GL1CrossingID;}
  short GetSpinGL1CrossingID() const {return SpinGL1CrossingID;}

  unsigned int GetTrigScaled() const {return TrigScaled;}
  unsigned int GetTrigLive()   const {return TrigLive;}


  void SetRunID(int rid)    {RunID = rid;}
  void SetEventID(int evid) {EventID = evid;}
  void SetZVertex(float vt) {VtxZ = vt;}

  void SetGL1CrossingID(short crossing)     {GL1CrossingID = crossing;}
  void SetSpinGL1CrossingID(short crossing) {SpinGL1CrossingID = crossing;}

  void SetTrigScaled(unsigned int trig) {TrigScaled = trig;}
  void SetTrigLive(unsigned int trig)   {TrigLive   = trig;}

  void Copy(IHeader *hdr);

protected:
  ClassDef(IHeader,1)
};
//=============================================================================

//inline functions


inline
void IHeader::Copy(IHeader *hdr)
{
  RunID   = hdr->RunID;
  EventID = hdr->EventID;
  VtxZ    = hdr->VtxZ;

  GL1CrossingID     = hdr->GL1CrossingID;
  SpinGL1CrossingID = hdr->SpinGL1CrossingID;

  TrigScaled = hdr->TrigScaled;
  TrigLive   = hdr->TrigLive;
}
//-----------------------------------------------------------------------------
#endif
