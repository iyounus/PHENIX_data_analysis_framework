
#ifndef IPI0RECO_HH
#define IPI0RECO_HH

#include <vector>

class IPdst;
class TH1D;

class IPi0Reco
{
public:
  IPi0Reco(int runNum);
  ~IPi0Reco();

private:
  int  _runNum;
  long _nEntries;

  IPdst *_pdst;


  std::vector <TH1D *> *_vHistos;

  TH1D **_hmass_PbSc;
  TH1D **_hpt_PbSc;
  TH1D **_hphi_PbSc;
  TH1D **_htheta_PbSc;

  TH1D **_hmass_PbGl;
  TH1D **_hpt_PbGl;
  TH1D **_hphi_PbGl;
  TH1D **_htheta_PbGl;

  TH1D **_hmass_ScGl;
  TH1D **_hpt_ScGl;
  TH1D **_hphi_ScGl;
  TH1D **_htheta_ScGl;

  short *_WarnMap;

  void CreateHistos();
  void Write();
  void ReadEMCalWarnMap(const char* file);
};
#endif
//-----------------------------------------------------------------------------
