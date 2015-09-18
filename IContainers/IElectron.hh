//-----------------------------------------------------------
// IElectron.hh
//	 Created  Feb 09, 2010
//       Imran Younus
//
//
// ***   IMPORTANT   ***
// This class stores extra variable required for electron analysis in central
// arm. These varialble are not stored in ITrack class because most of the
// time there is no RICH information. Only if n0 or n1 value is available
// the track is an electron candidate, and then this object is stored in the 
// tree. Therefore, this class is just an extenstion to ITrack.
// In order to associate this electron to the ITrack, varialbe Track is stored
// here. This is just an index, and tells you the track number. NOTE: counting
// starts from 0, so that it can simply be used as an index to an array. So,
// If there are 4 tracks in an event, this Track varialbe can be 0,1,2 or 3,
// depending on which track it is associated with. This is VERY IMPORTANT.
//-----------------------------------------------------------

#ifndef IELECTRON_H
#define IELECTRON_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif


class IElectron : public TObject
{
private:
  // See the description of Track variable above.
  short Track;        // associated track# (IMPORTANT: counting starts from 0).

  short N0;           // PHCentrakTrack::get_n0
  short N1;           // PHCentrakTrack::get_n1

  // This is the number of photo-electrons (i.e. summed pulse height) summed
  // over the normally sized ring area.
  float Npe0;         // PHCentralTrack::get_npe0
  float Chi2;         // PHCentralTrack::get_chi2

  // This variable measures the displacement of ring center w.r.t. the track
  // projection positon in the RICH PMT array.
  float RICHdisp;     // PHCentralTrack::get_disp


  // NOTE: The following EMCal matching variables are different than the ones
  // stored in ITrack. Here is the description of these from PHENIX web page
  // https://www.phenix.bnl.gov/WWW/offline/CNT.html
  // emcsdphi_e : The position resolution of the EMC depends upon the shower
  // type. This calculation determines the dphi in SIGMAS assuming the
  // resoution appropriate for EM showers.
  // emcsdz_e : The position resolution of the EMC depends upon the shower
  // type. This calculation determines the dz in SIGMAS assuming the resoution
  // appropriate for EM showers.

  float EmcsdPhiE;    // PHCentralTrack::get_emcsdphi_e
  float EmcsdZE;      // PHCentralTrack::get_emcsdz_e


public:
  IElectron(){;}
  ~IElectron(){;}


  short GetTrack() const {return Track;}

  short GetN0()    const {return N0;}
  short GetN1()    const {return N1;}

  float GetNpe0()  const {return Npe0;}
  float GetChi2()  const {return Chi2;}

  float GetRICHdisp()  const {return RICHdisp;}
  float GetEmcsdPhiE() const {return EmcsdPhiE;}
  float GetEmcsdZE()   const {return EmcsdZE;}


  void SetTrack(short tr) {Track = tr;}

  void SetN0(short n0)    {N0 = n0;}
  void SetN1(short n1)    {N1 = n1;}

  void SetNpe0(float npe) {Npe0 = npe;}
  void SetChi2(float chi) {Chi2 = chi;}

  void SetRICHdisp(float disp)   {RICHdisp = disp;}
  void SetEmcsdPhiE(float sdphi) {EmcsdPhiE = sdphi;}
  void SetEmcsdZE(float sdz)     {EmcsdZE = sdz;}

  void Copy(IElectron *electron);


protected:
  ClassDef(IElectron,1)
};
//=============================================================================

//inline functions


inline
void IElectron::Copy(IElectron *electron)
{
  Track = electron->Track;

  N0  = electron->N0;
  N1  = electron->N1;

  Npe0 = electron->Npe0;
  Chi2 = electron->Chi2;

  RICHdisp = electron->RICHdisp;
  EmcsdPhiE = electron->EmcsdPhiE;
  EmcsdZE = electron->EmcsdZE;
}
//-----------------------------------------------------------------------------
#endif
