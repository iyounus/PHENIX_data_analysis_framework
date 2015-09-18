
#ifndef PI0MASS_HH
#define PI0MASS_HH

class TH1D;
class TFile;
class TF1;
class TGraph;


class Pi0Mass
{
public:
  Pi0Mass(const char* instance, const char* det, const char* file);
  ~Pi0Mass();


  TH1D *Get_hmass(int i) {return _hmass[i];}
  TH1D *Get_hpt(int i)   {return _hpt[i];}

  TF1  *Get_sign(int i)  {return _sign[i];}
  TF1  *Get_bkgr(int i)  {return _bkgr[i];}

  TGraph *Get_gSB()      {return _gSB;}

  void Draw_hmass(int i);
  void Fit_hmass(int i, char *opt="R");
  TF1  *Fit_hMass(char *opt="R", double mean=0.135, double sig=0.9);
  TH1D *Get_hMass() {return _hMass;}

  double Get_nPi0s(int i){return _nPi0s[i];}
  double Get_SS(int i)   {return _SS[i];}
  double Get_BB(int i)   {return _BB[i];}
  double Get_SB(int i)   {return _SB[i];}

private:
  TFile *_f;

  TH1D *_hmass[8];
  TH1D *_hpt[8];
  TH1D *_hMass;
  TGraph *_gSB;

  TF1 *_fn[10];
  TF1 *_sign[10];
  TF1 *_bkgr[10];

  double _nPi0s[10];
  double _meanPt[10];
  double _SS[10];
  double _BB[10];
  double _SB[10];
};
#endif
