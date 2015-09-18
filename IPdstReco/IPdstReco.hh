//-----------------------------------------------------------
// IPdstReco.hh
//       Created April 11 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IPDSTRECO_HH
#define IPDSTRECO_HH

#include "SubsysReco.h"


class PHGlobal;
class EventHeader;
class PHCompositeNode;
class PHCentralTrack;
class PHMuoTracksOut;
class emcClusterContainer;
class emcClusterContent;
class ErtOut;
class SmdOut;
class SpinDataEventOut;
class TrigLvl1;
class lpcRaw;


class IHeaderList;
class IPhotonList;
class ITrackList;
class IElectronList;
class IMuonList;
class IDiMuonList;
class IBbcList;
class IZdcSmdList;


class TFile;
class TTree;


class IPdstReco : public SubsysReco
{
public:
  IPdstReco(char* outFile);
  virtual ~IPdstReco();

  int End(PHCompositeNode *topNode);
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode){ return 0; }
  const char *Name() const { return ThisName.c_str(); }
  void Print(const char *what = "ALL") const { return; }
  void Verbosity(const int ival) {verbosity = ival;}


  // These functions are used to set different cuts
  // If needed, these must be called before Fun4AllServer::registerSubsystem
  // The default values are also given here.
  void SetZVtxCentral    (float z = 35.)  {_cZVtxCentral = z;}
  void SetZVtxMuon       (float z = 55.)  {_cZVtxMuon = z;}
  void SetMaxZedDC       (float z = 80.)  {_cMaxZedDC = z;}
  void SetMinTrackPt     (float pt = 0.2) {_cMinTrackPt = pt;}
  void SetMinPhotonEn    (float en = 0.1) {_cMinPhotonEn = en;}
  void SetMinPhotonProb  (float p = 0.02) {_cMinPhotonProb = p;}
  void SetMinTrackHighPt (float pt = 0.2) {_cMinTrackHighPt = pt;}
  void SetMinPhotonHighEn(float en = 0.8) {_cMinPhotonHighEn = en;}


  // These functions toggle different flags to configure final output.
  // If needed, these must be called before Fun4AllServer::registerSubsystem
  void SetSpinFlag    (bool spin = true)  {_bSpin = spin;}
  void SetErtFlag     (bool ert  = true)  {_bErt = ert;}
  void SetMuTrigFlag  (bool mut  = false) {_bMuTrig = mut;}
  void SetMinBiasFlag (bool mb   = false) {_bMinBias = mb;}
  void SetPhotonFlag  (bool pho  = true)  {_bPhoton = pho;}
  void SetTrackFlag   (bool trk  = true)  {_bTrack = trk;}
  void SetMuonFlag    (bool muon = true)  {_bMuon = muon;}
  void SetDiMuonFlag  (bool muon = true)  {_bDiMuon = muon;}
  void SetBbcFlag     (bool bbc  = true)  {_bBbc = bbc;}
  void SetZdcSmdFlag  (bool zdc  = true)  {_bZdcSmd = zdc;}


protected:
  // data nodes to be read
  PHGlobal            *d_global;
  EventHeader         *d_evthdr;
  PHCentralTrack      *d_cnttrack;
  PHMuoTracksOut      *d_muotrack;
  emcClusterContainer *d_emcclust;
  ErtOut              *d_ertout;
  SmdOut              *d_smdout;
  SpinDataEventOut    *d_spindata;
  TrigLvl1            *d_triglvl1;
  lpcRaw              *d_lpcRaw;

private:
  char *_outFileName;
  TFile *_outFile;
  TTree *_pdst;

  //data nodes to write out
  IHeaderList   *d_header;
  IPhotonList   *d_photon;
  ITrackList    *d_track;
  IElectronList *d_elect;
  IMuonList     *d_muon;
  IDiMuonList   *d_dimuon;
  IBbcList      *d_bbc;
  IZdcSmdList   *d_zdcsmd;

  bool *_EMC_badChan;


  // cuts used in the pdsts
  float _cZVtxCentral;  // zvertex cut for centeral arm events
  float _cZVtxMuon;     // zvertex cut for muon arm events
  float _cMaxZedDC;     // max ZedDC cut
  float _cMinTrackPt;   // every track in an event should have pt>_MinTrackPt
  float _cMinPhotonEn;  // every photon in an event should have En>_MinPhotonEn

  float _cMinPhotonProb;   // photon probability cut
  float _cMinTrackHighPt;  // If there is at least one track in an event with
                           // pt > _MinTrackHighPt, all tracks are saved, 
                           // otherwise all the tracks are discarded.
  float _cMinPhotonHighEn; // If there is at least on photon in an event with
                           // En > _MinPhotonHighEn, all photons are saved,
                           // otherwise all the photons are discarded.



  // these flags are used to make code more configurable without commenting
  // out lines of code. By default all these are true.
  bool _bSpin;    // set false for heavy ion data set.
  bool _bErt;     // set false for MinBias, Muon or MPC trigger data set
  bool _bMuTrig;  // set true for muon trigger data
  bool _bMinBias; // set true for min bias trigger data

  bool _bPhoton;  // set false if no photons are needed.
  bool _bTrack;   // set false if no tracks are needed.
  bool _bMuon;    // set false if no muons are needed.
  bool _bDiMuon;  // set false if no di-muons are needed.
  bool _bBbc;     // set false if no bbc info is needed.
  bool _bZdcSmd;  // set false if no zdc info is needed.


  int _runNumber;
  int _eventNumber;
  int _nTotalEvts;

  //for verbosity
  int _ival;
  short _ErrCount;

  short FillTracks();
  short FillPhotons();
  short FillMuons();
  void  FillDiMuons(int n, int *imuon);

  void FillBbc();
  void FillZdcSmd();

  void InitTree();

  short TowerID(short arm, short sect, short ypos, short zpos);
  bool  CompareTowerID(short towerID);
  short GetEMCalSM(short arm, short sect, short ypos, short zpos);
  short GetERT(short arm, short sect, short ypos, short zpos);
  void  ReadEMCalWarnMap(const char* file);
};
//=============================================================================


inline
short IPdstReco::TowerID(short arm, short sect, short ypos, short zpos)
{
  short towerID=-99;

  if (ypos < 0 || zpos < 0 || sect < 0) return towerID;

  if (arm == 0)
    towerID = zpos + 72*ypos + sect*2592;

  if (arm == 1 && sect > 1)
    towerID = zpos + 72*ypos + (sect+2)*2592;

  if (arm == 1 && sect < 2)
    towerID = zpos + 96*ypos + sect*4608 + 15552;

  return towerID;
}
//-----------------------------------------------------------------------------


inline
short IPdstReco::GetEMCalSM(short arm, short sector, short ytwr, short ztwr)
{
  short sm;
  if(arm == 0 || sector >= 2)         // PbSc
    sm = (ytwr/12) * 6 + (ztwr/12);
  else
    sm = (ytwr/12) * 8 + (ztwr/12);   // PbGl
  
  return sm;
}
//-----------------------------------------------------------------------------

#endif
