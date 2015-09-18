
#ifndef IGAMMAHADRON_HH
#define IGAMMAHADRON_HH

#include <vector>

class TH1D;
class TFile;
class TTree;

class IPdst;
class IHeader;
class IPhoton;
class ITrack;

class IEvent;

class IGammaHadron
{
public:
  IGammaHadron();
  IGammaHadron(int run);
  ~IGammaHadron();

private:
  IPdst *_pdst;


  IEvent *_event;
  std::vector <IEvent *> *_vEvent;

  std::vector <TH1D *> *_vHist;


  int _nPtT;     // no. of ptt bins
  int _nPtA;     // no. of pta bins

  TH1D ***_hPtT; // 2D array of pointers. I cannot use IConsts to define array.
                 // thats why so many pointer.
  TH1D  **_hPtA; // associated pt is considered independent

  TH1D  **_hProb;// EMCal shower shape prob
  TH1D  **_hEP;  // energy/momentum  
  TH1D ***_hZt;

  TH1D   *_hPtT_all;
  TH1D   *_hPtA_all;

  TH1D  **_hEcone;
  TH1D  **_hDiPhoMs;  // di photon mass
  TH1D  **_hDiHadMs;  // di hadron mass
  TH1D ***_hDphiPosH;
  TH1D ***_hDphiNegH;


  TFile *_f;
  TTree *_gam;
  // These variables are used by the tree branches

  double _phoPhi;
  double _hadPhi;
  double _phoPt;
  double _hadPt;
  double _phoEn;
  double _EP;
  double _za;
  double _xE;
  double _dPhi;
  double _eCone;
  double _hadCh;
  double _diHadMs;
  double _diPhoMs;
  double _prob;



  void ReadFiles();
  void ReadFile(int run);

  void CreateTree();
  void CreateHistos();
  void FillHistos();
  void ResetVector();
};
#endif
//------------------------------------------------------------------------------
