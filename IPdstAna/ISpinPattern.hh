
#ifndef ISPINPATTERN_H
#define ISPINPATTERN_H

#include <map>

class TTree;
class TRandom3;
class TGraphErrors;
class ISpinPatternRun;


#define XINGs 120


class ISpinPattern
{
public:
  ISpinPattern(const char* spinInfoFile, const char* goodRunList,
	       const char* anaType="Run"); // "Run" or "Fill"
  ~ISpinPattern();

  enum SPIN {BPOS = +1, BNEG = -1, YPOS = +1, YNEG = -1,
	     BUP  = +1, BDN  = -1, YUP  = +1, YDN  = -1};

  enum RLUMI {R_ALL,           // L(++,--)/L(+-,-+)
	      R_PVIOLATING1,   // L(++)/L(--)
	      R_PVIOLATING2,   // L(+-)/L(-+)
	      R_ANBLUE,        // L(Bu)/L(Bd)
	      R1_ANBLUE,       // L(Bu,Yu)/L(Bu,Yd)
	      R2_ANBLUE,       // L(Bu,Yu)/L(Bd,Yu)
	      R3_ANBLUE,       // L(Bu,Yu)/L(Bd,Yd)
	      R_ANYELLOW,      // L(Yu)/L(Yd)
	      R1_ANYELLOW,     // L(Yu,Bu)/L(Yu,Bd)
	      R2_ANYELLOW,     // L(Yu,Bu)/L(Yd,Bu)
	      R3_ANYELLOW};    // L(Yu,Bu)/L(Yd,Bd)
  // NOTE: R1_BLUE = R2_YELLOW, R2_BLUE = R1_YELLOW, R3_BLUE = R3_YELLOW,


  void Shuffle();
  void KnuthShuffle();

  short GetSpinB(int run, int Xing);
  short GetSpinY(int run, int Xing);

  void GetSpinB(int run, short spin[XINGs]);
  void GetSpinY(int run, short spin[XINGs]);
  void GetBbcIN(int run, float bbc[XINGs]);

  short GetNBunches(int run, SPIN B, SPIN Y);
  short GetFillPattern(int run);

  double GetBbcSum(int run);
  double GetBbcSum(int run, SPIN B, SPIN Y);
  double GetBbcSum(SPIN B, SPIN Y);

  double GetBbcSum(int run, SPIN B1, SPIN Y1, SPIN B2, SPIN Y2);
  double GetBbcSum(SPIN B1, SPIN Y1, SPIN B2, SPIN Y2);

  void GetRelLumi(int run, RLUMI R_config, double &R, double &dR);
  void GetRelLumi(RLUMI R_config, double &R, double &dR);

  void GetPolarization(int run,
		       float &PolB, float &ePolBstat, float &ePolBsyst,
		       float &PolY, float &ePolYstat, float &ePolYsyst);

  void GetPolarizationB(int run, float &PolB,
			float &ePolBstat, float &ePolBsyst);

  void GetPolarizationY(int run, float &PolY,
			float &ePolYstat, float &ePolYsyst);



  TGraphErrors* RelLumiGraph(RLUMI R_config);

  bool Exists(int run);
  bool GoodRun(int run);

  void DrawBBC(int run);
  void DrawSpinB(int run);
  void DrawSpinY(int run);

  void SetSeed(int seed=0);
  void SetEvenBunches();    // this is to analyse even bunches only
  void SetOddBunches();     // this is to analyse odd  bunches only

  void Reset();

  void Print(int run);


private:
  int  _nGoodRuns;  // number of good runs
  int *_GoodRun;    // these are used when calculating bbc sum, rel lumi etc.

  TRandom3 *_RndmGen;
  std::map <int, ISpinPatternRun*> _SpinPat;

  bool _Rndm;
  bool _Even;
  bool _Odd;
  bool _FillAna;

  void ReadGoodRunsList(const char* file);
  void SetSpin_RunByRun(TTree *);
  void SetSpin_FillByFill(TTree *);
};
//----------------------------------------------------------------------------


class ISpinPatternRun
{
public:
  ISpinPatternRun(int runNum, short *spinB, short *spinY,
		  float *bbc_in, float *zdc_in);
  ~ISpinPatternRun();

  void Reset();
  void Shuffle(TRandom3 *rnd);
  void KnuthShuffle(TRandom3 *rnd);

  int   GetRunNum(){return _RunNum;}

  short GetFillPattern(){return _FillPat;}

  short GetSpinB(int Xing) {return _SpinB[Xing];}
  short GetSpinY(int Xing) {return _SpinY[Xing];}

  float GetBbcIn(int Xing);
  float GetZdcIn(int Xing);

  double GetBbcSum();
  double GetBbcSum(short B, short Y);
  double GetBbcSum(short B1, short Y1, short B2, short Y2);

  short GetNBunches(short B, short Y);

  void SetPolarization(float polB, float ePolBstat, float ePolBsyst,
		       float polY, float ePolYstat, float ePolYsyst);

  void GetPolarization(float &polB, float &ePolBstat, float &ePolBsyst,
		       float &polY, float &ePolYstat, float &ePolYsyst);

  void GetPolarizationB(float &polB, float &ePolBstat, float &ePolBsyst);
  void GetPolarizationY(float &polY, float &ePolYstat, float &ePolYsyst);



  void SetFillPattern(short *spinB, short *spinY);

  void SetEvenBunches();  // this is to analyse even bunches only
  void SetOddBunches();   // this is to analyse odd bunches only

  void Print();

private:
  int    _RunNum; // run number or fill number.

  float  _PolB;
  float  _PolY;

  float  _ePolBstat;   // these error are for a fill and not for a run
  float  _ePolBsyst;
  float  _ePolYstat;
  float  _ePolYsyst;

  short *_SpinBdb;     // these are actual spin values read from the file, and
  short *_SpinYdb;     // are not chagned at all.

  short *_SpinB;       // these values are used for different type of analyses
  short *_SpinY;       // e.g., even/odd bunches, random spins etc.
                       // bbc sum and rel lumi are calculated using these 
                       // spin signs.
  float *_BbcIn;
  float *_ZdcIn;

  short _FillPat;
};
#endif
//----------------------------------------------------------------------------
