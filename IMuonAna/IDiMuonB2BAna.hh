#ifndef IDIMUONB2BANA_HH
#define IDIMUONB2BANA_HH

#include <vector>

class IPdst;
class IEvent;
class IFillLookup;
class ISpinPattern;

class TH1D;
class TFile;
class TTree;


class IDiMuonB2BAna
{
 public:
  IDiMuonB2BAna(const char* outFile);
  ~IDiMuonB2BAna();

  enum SPIN  {Up=+1, Dn=-1};
  enum DET   {North, South};
  enum SIDE  {East, West};
  enum BEAM  {Blue, Yelo};
  enum SIGN  {MM, PM, PP};  // minus-minus, plus-minus, plus-plus

 private:
  IPdst *_pdst;
  TFile *_outFile;
  TTree *_dimu;

  std::vector<IEvent *> *_vEvent;
  std::vector<TH1D *>   *_vHist;

  ISpinPattern *_spin;
  short _spinB[120];
  short _spinY[120];

  IFillLookup *_fillLookup;

  TH1D *_hDiMuMass[3]; // for three differen charge sign combinations
  TH1D *_hDiMuPt[3];


  // tree variables
  int    _RunNumber;
  int    _FillNumber;
  int    _EvtNumber;
  short  _NMuons ;
  short  _SpinXingID;
  short  _SpinB;
  short  _SpinY;
  double _Zvtx;

  double _Mu1_DDG0;
  double _Mu1_DG0;
  double _Mu1_DS3;
  double _Mu1_DS3ctp;
  double _Mu1_MuTrChi2;
  double _Mu1_MuIdChi2;
  double _Mu1_Px;
  double _Mu1_Py;
  double _Mu1_Pz;
  double _Mu1_Pt;
  double _Mu1_P;
  short  _Mu1_Charge;
  short  _Mu1_nMuTrHits;
  short  _Mu1_nMuIdHits;
  unsigned int _Mu1_MuTrHitPat;
  unsigned int _Mu1_MuIdHitPat;

  double _Mu2_DDG0;
  double _Mu2_DG0;
  double _Mu2_DS3;
  double _Mu2_DS3ctp;
  double _Mu2_MuTrChi2;
  double _Mu2_MuIdChi2;
  double _Mu2_Px;
  double _Mu2_Py;
  double _Mu2_Pz;
  double _Mu2_Pt;
  double _Mu2_P;
  short  _Mu2_Charge;
  short  _Mu2_nMuTrHits;
  short  _Mu2_nMuIdHits;
  unsigned int _Mu2_MuTrHitPat;
  unsigned int _Mu2_MuIdHitPat;

  short  _diMuCharge;
  double _diMuMass;
  double _diMuP;
  double _diMuPt;
  double _diMuPz;
  double _diMuPhi;
  double _diMuEta;
  // --------------

  short  _Arm;
  short  _Side;
  short  _PtBin;

  void CreateHistos();
  void ReadPdst(int run);
  void FillHistos();
  void ResetEventVector();
};
#endif
