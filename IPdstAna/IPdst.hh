
#ifndef IPDST_HH
#define IPDST_HH


class TChain;

class IHeader;
class IBbc;
class IZdcSmd;
class IPhoton;
class ITrack;
class IElectron;
class IMuon;
class IDiMuon;

class IHeaderList;
class IBbcList;
class IZdcSmdList;
class IPhotonList;
class ITrackList;
class IElectronList;
class IMuonList;
class IDiMuonList;


class IPdst
{
public:
  // NOTE: TConsts must be defined before instantiating this class
  // IPdst requires IConsts::InDir and IConsts::FilePrefix

  IPdst(int run);             // this ctor takes only one run
  IPdst(short n, int *run);   // this ctor takes an array of run numbers,
                              // this can be used to combine runs in a fill.
  IPdst(const char* runList); // this ctor takes text file with list or runs 


  ~IPdst();

  void GetEntry(int i);
  int  GetEntries();
  void SetBranchStatus(const char *branchname, const int status);
  int  GetRunNumberFromPDST();

  bool BbcExists()    {return bool(_bbc)    ;}
  bool ZdcSmdExists() {return bool(_zdcsmd) ;}

  unsigned short GetNPhotons()   {return _nPhotons  ;}
  unsigned short GetNTracks()    {return _nTracks   ;}
  unsigned short GetNElectrons() {return _nElectrons;}
  unsigned short GetNMuons()     {return _nMuons    ;}
  unsigned short GetNDiMuons()   {return _nDiMuons  ;}
  unsigned short GetNNeutrons()  {return _nNeutrons ;}


  IHeader* GetHeader()         {return _header;}
  IBbc*    GetBbc()            {return _bbc   ;}
  IZdcSmd* GetZdcSmd()         {return _zdcsmd;}

  IPhoton*   GetPhoton  (short i) {return _photon[i];}
  ITrack*    GetTrack   (short i) {return _track[i] ;}
  IElectron* GetElectron(short i) {return _elect[i] ;}
  IMuon*     GetMuon    (short i) {return _muon[i]  ;}
  IDiMuon*   GetDiMuon  (short i) {return _dimuon[i];}


  double DiPhotonMass(short i, short j);// di photon mass for ith and jth photon
  double DiHadronMass(short i, short j);// di hadron mass for ith and jth track
  double DiMuonMass  (short i, short j);// di muon   mass for ith and jth muon

  // energy in a cone of ConeAngle around ith photon
  double ECone(unsigned short i, double ConeAngle=0.5);

  bool CompareTowerID(short tw);

private:
  TChain *_pDST;

  unsigned short _nPhotons;
  unsigned short _nTracks;
  unsigned short _nElectrons;
  unsigned short _nMuons;
  unsigned short _nDiMuons;
  unsigned short _nNeutrons;

  IHeader    *_header;
  IBbc       *_bbc;
  IZdcSmd    *_zdcsmd;
  IPhoton   **_photon;
  ITrack    **_track;
  IElectron **_elect;
  IMuon     **_muon;
  IDiMuon   **_dimuon;


  IHeaderList   *_headerList;
  IBbcList      *_bbcList;
  IZdcSmdList   *_zdcsmdList;
  IPhotonList   *_photonList;
  ITrackList    *_trackList;
  IElectronList *_electList;
  IMuonList     *_muonList;
  IDiMuonList   *_dimuonList;

  void SetBranches();
};
#endif
