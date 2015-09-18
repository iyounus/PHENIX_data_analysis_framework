
#ifndef IPDSTSHORT_HH
#define IPDSTSHORT_HH

class TFile;
class TTree;

class IPdst;
class IHeaderList;
class IPhotonList;
class ITrackList;
class IElectronList;
class IMuonList;
class IBbcList;
class IZdcSmdList;

typedef unsigned short ushort;


class IPdstShort
{
public:
  IPdstShort(int run, const char* prefix = "IPdstShort_");
  ~IPdstShort();

  void Process();


  // These functions are used to set different cuts
  // The default values are also given here.
  void SetMinNPhotons    (ushort n = 1)   {_cMinNPhotons = n;}
  void SetMinNTracks     (ushort n = 0)   {_cMinNTracks = n;}
  void SetMinNMuons      (ushort n = 0)   {_cMinNMuons = n;}
  void SetZVtxCentral    (float z = 35.)  {_cZVtxCentral = z;}
  void SetZVtxMuon       (float z = 55.)  {_cZVtxMuon = z;}
  void SetMaxZedDC       (float z = 75.)  {_cMaxZedDC = z;}
  void SetMinZedDC       (float z = 5.0)  {_cMaxZedDC = z;}
  void SetMinTrackPt     (float pt = 0.2) {_cMinTrackPt = pt;}
  void SetMaxTrackPt     (float pt = 50.) {_cMaxTrackPt = pt;}
  void SetMinPhotonEn    (float en = 0.1) {_cMinPhotonEn = en;}
  void SetMaxPhotonEn    (float en = 50.) {_cMaxPhotonEn = en;}
  void SetMinPhotonProb  (float p = 0.02) {_cMinPhotonProb = p;}
  void SetMinTrackHighPt (float pt = 0.2) {_cMinTrackHighPt = pt;}
  void SetMinPhotonHighPt(float pt = 0.4) {_cMinPhotonHighPt = pt;}
  void SetMinEmcTof      (float tof =-10.){_cMinEmcTof = tof;}
  void SetMaxEmcTof      (float tof = 15.){_cMaxEmcTof = tof;}

  void SetPhotonFlag  (bool pho = true) {_bPhoton = pho;}
  void SetTrackFlag   (bool trk = true) {_bTrack = trk;}
  void SetElectronFlag(bool ele = true) {_bElect = ele;}
  void SetMuonFlag    (bool muo = true) {_bMuon = muo;}
  void SetBbcFlag     (bool bbc = true) {_bBbc = bbc;}
  void SetZdcSmdFlag  (bool zdc = true) {_bZdcSmd = zdc;}

  void SetEMCalWarnMap(const char* warnmap);

private:
  TFile *_NewFile;
  TTree *_NewPdst;

  IPdst *_pdst;


  // cuts used in the pdsts
  ushort _cMinNPhotons;  // min number of photons required in event.
  ushort _cMinNTracks;   // min number of tracks required in event.
  ushort _cMinNMuons;    // min number of muons required in event.

  float _cZVtxCentral;  // zvertex cut for centeral arm events
  float _cZVtxMuon;     // zvertex cut for muon arm events
  float _cMaxZedDC;     // max ZedDC cut
  float _cMinZedDC;     // min ZedDC cut
  float _cMinTrackPt;   // every track in an event should have Pt>_MinTrackPt
  float _cMaxTrackPt;   // every track in an event should have Pt<_MaxTrackPt
  float _cMinPhotonEn;  // every photon in an event should have Pt>_MinPhotonEn
  float _cMaxPhotonEn;  // every photon in an event should have Pt<_MaxPhotonEn

  float _cMinPhotonProb;   // photon probability cut
  float _cMinTrackHighPt;  // If there is at least one track in an event with
                           // Pt > _MinTrackHighPt, all tracks are saved, 
                           // otherwise all the tracks are discarded.
  float _cMinPhotonHighPt; // If there is at least on photon in an event with
                           // Pt > _MinPhotonHighEn, all photons are saved,
                           // otherwise all the photons are discarded.
  float _cMinEmcTof;
  float _cMaxEmcTof;


  // these flags are used to turn on/off different branches in output pDST
  bool _bPhoton;  // set false if no photons are needed.
  bool _bTrack;   // set false if no tracks are needed.
  bool _bElect;   // set false if no electrons are needed.
  bool _bMuon;    // set false if no muons are needed.
  bool _bBbc;     // set false if no bbc info is needed.
  bool _bZdcSmd;  // set false if no zdc info is needed.


  bool _ApplyWarnMap; // if true, then EMCal warn map will be applied
  bool *_EMC_badChan; // list of bad channels.


  // these are the objects that can be written to the output file if
  // the branch is created.
  IHeaderList   *_lHeader;
  IPhotonList   *_lPhoton;
  ITrackList    *_lTrack;
  IElectronList *_lElect;
  IMuonList     *_lMuon;
  IBbcList      *_lBbc;
  IZdcSmdList   *_lZdcSmd;

  void Init();

  ushort FillPhotons();
  ushort FillTracks();
  ushort FillMuons();
};
#endif
