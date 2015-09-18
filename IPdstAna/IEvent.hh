#ifndef IEVENT_HH
#define IEVENT_HH

#include <vector>

class IHeader;
class IBbc;
class IZdcSmd;
class IPhoton;
class ITrack;
class IElectron;
class IMuon;
class IDiMuon;


class IEvent
{
public:
  IEvent();
  ~IEvent();

  unsigned short GetNPhotons()   {return _photon->size();}
  unsigned short GetNTracks()    {return _track->size() ;}
  unsigned short GetNElectrons() {return _elect->size() ;}
  unsigned short GetNMuons()     {return _muon->size()  ;}
  unsigned short GetNDiMuons()   {return _dimu->size()  ;}

  unsigned short GetNNeutrons()  {return bool(_zdcsmd) ? 1 : 0;}


  IHeader*   GetHeader() {return _header;}
  IBbc*      GetBbc()    {return _bbc   ;}
  IZdcSmd*   GetZdcSmd() {return _zdcsmd;}

  IPhoton*   GetPhoton  (short i) {return (*_photon)[i];}
  ITrack*    GetTrack   (short i) {return (*_track)[i] ;}
  IElectron* GetElectron(short i) {return (*_elect)[i] ;}
  IMuon*     GetMuon    (short i) {return (*_muon)[i]  ;}
  IDiMuon*   GetDiMuon  (short i) {return (*_dimu)[i]  ;}


  void AddHeader(IHeader *hdr);
  void AddBbc   (IBbc    *bbc);
  void AddZdcSmd(IZdcSmd *zdc);

  void AddPhoton  (IPhoton   *pho);
  void AddTrack   (ITrack    *trk);
  void AddElectron(IElectron *ele);
  void AddMuon    (IMuon     *muo);
  void AddDiMuon  (IDiMuon   *dmu);


  double DiPhotonAngle(short i, short j);// angle between ith and jth photons
  double DiPhotonMass (short i, short j);// inv-mass for ith and jth photon
  double DiHadronMass (short i, short j);// inv-mass for ith and jth track
  double DiMuonMass   (short i, short j);// inv-mass for ith and jth muon

  // energy in a cone of ConeAngle around ith photon
  double ECone(unsigned short i, double ConeAngle=0.5);


private:
  IHeader  *_header;
  IBbc     *_bbc;
  IZdcSmd  *_zdcsmd;

  std::vector<IPhoton*>   *_photon;
  std::vector<ITrack*>    *_track;
  std::vector<IElectron*> *_elect;
  std::vector<IMuon*>     *_muon;
  std::vector<IDiMuon*>   *_dimu;
};
#endif
