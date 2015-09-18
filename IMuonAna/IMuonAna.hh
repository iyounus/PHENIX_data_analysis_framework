#ifndef IMUONANA_HH
#define IMUONANA_HH

#include <vector>

class IPdst;
class IEvent;
class TH1D;
class TFile;

class IMuonAna
{
public:
  IMuonAna();
  ~IMuonAna();


private:
  TFile *_outFile;

  IPdst *_pdst;

  std::vector <IEvent *> *_vEvent;
  std::vector <TH1D *> *_vHist;

  enum DET {N, S};

  TH1D *_hPt[2][4];
  TH1D *_hPhi[2][4];
  TH1D *_hPz[2][4];
  TH1D *_hMuIdHits;
  TH1D *_hMuTrHits;
  TH1D *_hDG0;
  TH1D *_hDDG0;
  TH1D *_hDS3;
  TH1D *_hDS3cpt;

  void CreateHistos();
  void ReadPdst(int run);
  void FillHistos();
  void ResetEventVector();
  void ResetHistos();
};
#endif
