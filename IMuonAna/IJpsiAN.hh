#ifndef IJPSIAN_HH
#define IJPSIAN_HH


class TF1;
class TH1D;
class TFile;
class TGraph;

class ISpinPattern;

#define nSPINCONF   16
#define nASYMCONF1  8
#define nASYMCONF2  14
#define nBEAM       2  // blue and yellow
#define nARM        2  // north and south
#define nSIDE       2  // east and west

 

class IJpsiAN
{
public:
  IJpsiAN();
  ~IJpsiAN(){;}


  // Here left means y x p > 0, where y is the unit vector along
  // y-axis in PHENIX lab frame, and p is the momentum direction of
  // the beam. Similarly right is y x p < 0. 

  enum SPINCONF {Blue_Up_Left_Forward,   Blue_Dn_Left_Forward,
		 Blue_Up_Right_Forward,  Blue_Dn_Right_Forward,
		 Blue_Up_Left_Backward,  Blue_Dn_Left_Backward,
		 Blue_Up_Right_Backward, Blue_Dn_Right_Backward,
		 Yelo_Up_Left_Forward,   Yelo_Dn_Left_Forward,
		 Yelo_Up_Right_Forward,  Yelo_Dn_Right_Forward,
		 Yelo_Up_Left_Backward,  Yelo_Dn_Left_Backward,
		 Yelo_Up_Right_Backward, Yelo_Dn_Right_Backward};


  // asymmetry if (Up - Dn)/(Up + Dn).
  enum ASYMCONF {Blue_Left_Forward,  Blue_Right_Forward,
		 Blue_Left_Backward, Blue_Right_Backward,
		 Yelo_Left_Forward,  Yelo_Right_Forward,
		 Yelo_Left_Backward, Yelo_Right_Backward,
		 Blue_Forward, Blue_Backward,
		 Yelo_Forward, Yelo_Backward,
		 Forward, Backward};

  enum SP    {Up = 1, Dn = -1};
  enum DET   {North, South};
  enum SIDE  {East, West};
  enum BEAM  {Blue, Yelo};


  enum PAR {A, K, Ajpsi, Mjpsi, Sjpsi, Apsi, Mpsi, Spsi};


private:
  TFile *_file;

  ISpinPattern *_spin;

  double *_Fill;
  int  _nFills;

  int _NPtBins;

  double _JpsiMs[nARM];
  double _JpsiSig[nARM];
  double _PsiMs[nARM];

  double ***_Njpsi;    // array [PtBins][nSPINCONF][FillNum]
  double ***_Nbkgr;    // array [PtBins][nSPINCONF][FillNum]

  double ***_AN;       // array [PtBins][nASYMCONF2][FillNum]
  double ***_dAN;      // array [PtBins][nASYMCONF2][FillNum]

  double ***_SumPhi;
  double ***_SumCosPhi;
  double ***_AccFunc;
  double ***_dAccFunc;


  float  *_Pol[nBEAM];
  float  *_ePolStat[nBEAM];
  float  *_ePolSyst[nBEAM];
  double *_R[nBEAM];
  double *_dR[nBEAM];


  TH1D *_hJpsiMs[nARM];
  TH1D *_hDiMuMs[nARM];

  TF1 *_fMs[nARM];

  TGraph *(*_gNjpsi)[nSPINCONF];

  void CreateArrays();
  void GetSpinInfo();
  void GetHistos();
  void SetMassFunc();
  void CalculateAccFunc();
  void CalculateAsym();
  void CountJpsi();
  int  BlueSpinConfig(double phi, short arm, short spin);
  int  YeloSpinConfig(double phi, short arm, short spin);
};
#endif
