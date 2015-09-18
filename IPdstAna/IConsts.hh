
#ifndef ICONSTS_HH
#define ICONSTS_HH


#define  EPS          1.0e-20
#define  PI           3.14159265358979323846
#define  TwoPI        6.28318530717958647692
#define  PIby2        1.57079632679489661923
#define  SqrtTwoPI    2.50662827463100050241



class IConsts
{
public:
  IConsts(int run);
  ~IConsts();


  static int     RUN;

  static int     NPtTBins;   // number of trigger particle pt bins
  static int     NPtABins;   // number of associated particle pt bins
  static double* PtT;        // trigger particle pt bins
  static double* PtA;        // associated particle pt bins

  static double  MinPhi;
  static double  MaxPhi;
  static double  MinDPhi;
  static double  MaxDPhi;
  static double  MinDPhiN;   // min value of near side dphi histogram 
  static double  MaxDPhiN;   // max value of near side dphi histogram 
  static double  MinDPhiF;   // min value of away side dphi histogram 
  static double  MaxDPhiF;   // max value of away side dphi histogram 
  static double  MinTheta;
  static double  MaxTheta;

  static double  EMCB1;      // min phi for EMCal
  static double  EMCB2;      // max phi for EMCal
  static double  DCHB1;      // min phi for DCH
  static double  DCHB2;      // max phi for DCH

  static double  NSigN;      // no. of sigmas for near side dphi histo limits
  static double  NSigF;      // no. of sigmas for away side dphi histo limits
  static double* SigN;       // sigma for near side peak
  static double* SigF;       // sigma for away side peak

  static double  EAsym;
  static double  Pi0Width;   // sigma for pi0 mass width
  static double  ConeAngle;  // angle of cone arong a particle
  static double  MatchSig;   // PC3 & EMC matching sigma

  static short   TrChar;     // track charge. if 0 both + and - are included
  static short   EMCArm;     // +1 for west, -1 for east, 0 for both arms

  static double  BbcChLo;    // BBC charge lo limit
  static double  BbcChHi;    // BBC charge hi limit
  static short   BbcNLo;     // BBC multiplicity lo limit
  static short   BbcNHi;     // BBC multiplicity hi limit

  static short   AsymConfig; // asymmetry configuration:
                             // 0 -> physics asym, (++,--)-(+-,-+)
                             // 1 -> parity violating, (++)-(--)
                             // 2 -> parity violating, (+-)-(-+)
                             // 3 -> parity violating, (B+)-(B-)
                             // 4 -> parity violating, (Y+)-(Y-)
  static short   H[10][4];   // helicity combinations. see implementation

  static char*   RunList;    // txt file with list of runs to be processed
  static char*   InDir;      // directory with input pdst files
  static char*   OutDir;     // directory for any output files
  static char*   OutFile;    // output file name if any
  static char*   SpinFile;   // root file contining spin information
  static char*   PbScStats;  // root file with pi0 mass/widths for PbSc
  static char*   PbGlStats;  // root file with pi0 mass/widths for PbGl
  static char*   EMCWarnMap; // warn/dead map file for EMCal
  static char*   FillTable;  // fill lookup take for runs
  static char*   FilePrefix; // prefix for pdst file name: IPdstReco_


  static bool    Defined;    // IConsts is instantiated if true
  //void Read(const char* txt_file);
  static void Print();


  void SetRunList   (char* list);
  void SetInDir     (char* dir );
  void SetOutDir    (char* dir );
  void SetOutFile   (char* file);
  void SetSpinFile  (char* spin);
  void SetPbScStats (char* pbsc);
  void SetPbGlStats (char* pbsc);
  void SetEMCWarnMap(char* map );
  void SetFillTable (char* fill);
  void SetFilePrefix(char* pre );


  void SetTriggerPtBins(int nbins, double *ptt);
  void SetAssociatedPtBins(int nbins, double *pta);

  void SetAsymConfig(short a){AsymConfig = a;}
  void SetTrCharge(short ch){TrChar = ch;}
  void SetEMCalArm(short arm){EMCArm = arm;}
  void SetBBCchar(double lo, double hi){BbcChLo = lo; BbcChHi = hi;}
  void SetBBCmult(short  lo, short  hi){BbcNLo  = lo; BbcNHi  = hi;}
  void SetPhiAcceptanceT(double lo, double hi)
  {EMCB1 = lo*PI/180.; EMCB2 = hi*PI/180.;}
  void SetPhiAcceptanceA(double lo, double hi)
  {DCHB1 = lo*PI/180.; DCHB2 = hi*PI/180.;}
  void SetPi0Width(double Nsigma){Pi0Width = Nsigma;}
  void SetPC3EMCmatch(double sig){MatchSig = sig;}

  void SetConeAngle(double angle){ConeAngle = angle;}

};
#endif
