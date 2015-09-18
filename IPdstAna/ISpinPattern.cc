#include "ISpinPattern.hh"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TString.h"
#include "TRandom3.h"
#include "TGraphErrors.h"


using namespace std;

typedef map <int, ISpinPatternRun*> MapSpinPat;

ISpinPattern::ISpinPattern(const char* spinInfoFile, const char* goodRunList,
			   const char* type):
  _Rndm(false),
  _Even(false),
  _Odd(false)
{
  cout << "ISpinPattern::ISpinPattern" << endl;

  _RndmGen = new TRandom3(0);

  _GoodRun = new int[2000]; // I don't expect more that 2000 good runs.
  for (int i=0; i<2000; i++) _GoodRun[i] = -99;
  ReadGoodRunsList(goodRunList);

  TFile *f = new TFile(spinInfoFile);
  TTree *t = (TTree *)f->Get("T");

  if (strcmp(type, "Run") == 0)
    SetSpin_RunByRun(t);
  else if (strcmp(type, "Fill") == 0)
    SetSpin_FillByFill(t);
  else
    {
      cout << "Incorrect anaysis type " << type
	   << "\nPossible valuse are:\n"
	   << "\"Run\"   for run by run analysis\n"
	   << "\"Fill\"  for fill by fill analysis"
	   << endl;
      exit(0);
    }

  f->Close();
  delete f;

  //cout << _SpinPat.size() << endl;
}
//-----------------------------------------------------------------------------


ISpinPattern::~ISpinPattern()
{
  MapSpinPat::const_iterator it;

  for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
    delete it->second;
}
//-----------------------------------------------------------------------------


void ISpinPattern::SetSpin_RunByRun(TTree *t)
{
  cout << "ISpinPattern::SetSpin_RunByRun" << endl;

  int   runNum;
  int   filNum;

  float  PolB;
  float  ePolBstat;
  float  ePolBsyst;
  float  PolY;
  float  ePolYstat;
  float  ePolYsyst;

  short spinB[XINGs];
  short spinY[XINGs];
  float bbc_in[XINGs];
  float zdc_in[XINGs];

  t->SetBranchAddress("RunNumber",  &runNum);
  t->SetBranchAddress("FillNumber", &filNum);

  t->SetBranchAddress("PolB",       &PolB);
  t->SetBranchAddress("ePolBstat",  &ePolBstat);
  t->SetBranchAddress("ePolBsyst",  &ePolBsyst);
  t->SetBranchAddress("PolY",       &PolY);
  t->SetBranchAddress("ePolYstat",  &ePolYstat);
  t->SetBranchAddress("ePolYsyst",  &ePolYsyst);

  t->SetBranchAddress("SpinB",      spinB);
  t->SetBranchAddress("SpinY",      spinY);
  t->SetBranchAddress("BBCin",      bbc_in);
  t->SetBranchAddress("ZDCin",      zdc_in);


  for (int i=0; i<t->GetEntries(); i++)
    {
      t->GetEntry(i);

      if (!GoodRun(runNum)) continue;

      ISpinPatternRun *pat = 
	new ISpinPatternRun(runNum, spinB, spinY, bbc_in, zdc_in);

      pat->SetPolarization(PolB, ePolBstat, ePolBsyst,
			   PolY, ePolYstat, ePolYsyst);

      _SpinPat.insert(MapSpinPat::value_type(runNum, pat));
    }
}
//-----------------------------------------------------------------------------


void ISpinPattern::SetSpin_FillByFill(TTree *t)
{
  cout << "ISpinPattern::SetSpin_FillByFill" << endl;

  int   runNum;
  int   filNum;

  float  PolB = 0.;
  float  ePolBstat = 0.;
  float  ePolBsyst = 0.;
  float  PolY = 0.;
  float  ePolYstat = 0.;
  float  ePolYsyst = 0.;

  float  PolB1 = 0.;
  float  ePolBstat1 = 0.;
  float  ePolBsyst1 = 0.;
  float  PolY1 = 0.;
  float  ePolYstat1 = 0.;
  float  ePolYsyst1 = 0.;

  short spinB[XINGs];
  short spinY[XINGs];
  float bbc_in[XINGs];
  float zdc_in[XINGs];

  short spinB_fill[XINGs];
  short spinY_fill[XINGs];
  float bbc_in_fill[XINGs];
  float zdc_in_fill[XINGs];

  for (int i=0; i<XINGs; i++)
    {
      spinB_fill[i]  = 0;
      spinY_fill[i]  = 0;
      bbc_in_fill[i] = 0;
      zdc_in_fill[i] = 0;
    }

  t->SetBranchAddress("RunNumber",  &runNum);
  t->SetBranchAddress("FillNumber", &filNum);

  t->SetBranchAddress("PolB",       &PolB);
  t->SetBranchAddress("ePolBstat",  &ePolBstat);
  t->SetBranchAddress("ePolBsyst",  &ePolBsyst);
  t->SetBranchAddress("PolY",       &PolY);
  t->SetBranchAddress("ePolYstat",  &ePolYstat);
  t->SetBranchAddress("ePolYsyst",  &ePolYsyst);

  t->SetBranchAddress("SpinB",      spinB);
  t->SetBranchAddress("SpinY",      spinY);
  t->SetBranchAddress("BBCin",      bbc_in);
  t->SetBranchAddress("ZDCin",      zdc_in);

  int oldFill = 0;
  for (int i=0; i<t->GetEntries(); i++)
    {
      t->GetEntry(i);
      if (!GoodRun(runNum)) continue;

      if (oldFill == 0) oldFill = filNum;

      if (filNum != oldFill)
	{
	  ISpinPatternRun *pat =
	    new ISpinPatternRun(oldFill, spinB_fill, spinY_fill,
				bbc_in_fill, zdc_in_fill);

	  pat->SetPolarization(PolB1, ePolBstat1, ePolBsyst1,
			       PolY1, ePolYstat1, ePolYsyst1);

	  _SpinPat.insert(MapSpinPat::value_type(oldFill, pat));

	  for (int i=0; i<XINGs; i++)
	    {
	      bbc_in_fill[i] = 0;
	      zdc_in_fill[i] = 0;
	    }

	  oldFill = filNum;
	}


      PolB1      = PolB;
      ePolBstat1 = ePolBstat;
      ePolBsyst1 = ePolBsyst;
      PolY1      = PolY;
      ePolYstat1 = ePolYstat;
      ePolYsyst1 = ePolYsyst;

      for (int j=0; j<XINGs; j++)
	{
	  spinB_fill[j] = spinB[j];
	  spinY_fill[j] = spinY[j];

	  bbc_in_fill[j] += bbc_in[j];
	  zdc_in_fill[j] += zdc_in[j];
	}
    }

  ISpinPatternRun *pat =
    new ISpinPatternRun(oldFill, spinB_fill, spinY_fill,
			bbc_in_fill, zdc_in_fill);

  pat->SetPolarization(PolB, ePolBstat, ePolBsyst,
		       PolY, ePolYstat, ePolYsyst);

  _SpinPat.insert(MapSpinPat::value_type(oldFill, pat));


//   cout << "In spin_patter.cc ------------------------------------" << endl;
//   cout << _SpinPat.size() << endl;
//   MapSpinPat::const_iterator it;
//   for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
//     cout << it->first << "\t" << it->second->_RunNum << endl;
//   cout << "In spin_patter.cc ------------------------------------" << endl;
}
//-----------------------------------------------------------------------------


bool ISpinPattern::Exists(int run)
{
  MapSpinPat::const_iterator it;
  it = _SpinPat.find(run);

  if (it == _SpinPat.end())
    {
      //cout << run << " not found" << endl;
      return false;
    }
  else
    return true;
}
//-----------------------------------------------------------------------------


void ISpinPattern::SetSeed(int seed)
{
  _RndmGen->SetSeed(seed);
}
//-----------------------------------------------------------------------------


void ISpinPattern::Shuffle()
{
  _RndmGen->SetSeed(0);

  MapSpinPat::const_iterator it;

  for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
    it->second->Shuffle(_RndmGen);

  _Rndm = true;
}
//-----------------------------------------------------------------------------


void ISpinPattern::KnuthShuffle()
{
  _RndmGen->SetSeed(0);

  MapSpinPat::const_iterator it;

  for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
    it->second->KnuthShuffle(_RndmGen);

  _Rndm = true;
}
//-----------------------------------------------------------------------------


void ISpinPattern::SetEvenBunches()
{
  cout << "ISpinPattern::SetEvenBunches" << endl;
  cout << "*********************************************************" << endl;
  cout << "***   Analysing EVEN bunches ONLY, don't forget it!   ***" << endl;
  cout << "*********************************************************" << endl;

  MapSpinPat::const_iterator it;

  for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
    it->second->SetEvenBunches();

  _Even = true;
}
//-----------------------------------------------------------------------------


void ISpinPattern::SetOddBunches()
{
  cout << "ISpinPattern::SetOddBunches" << endl;
  cout << "********************************************************" << endl;
  cout << "***   Analysing ODD bunches ONLY, don't forget it!   ***" << endl;
  cout << "********************************************************" << endl;

  MapSpinPat::const_iterator it;

  for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
    it->second->SetOddBunches();

  _Odd = true;
}
//-----------------------------------------------------------------------------


void ISpinPattern::Reset()
{
  MapSpinPat::const_iterator it;

  for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
    it->second->Reset();

  _Rndm = false;
  _Even = false;
  _Odd  = false;
}
//-----------------------------------------------------------------------------


short ISpinPattern::GetFillPattern(int run)
{
  return _SpinPat[run]->GetFillPattern();
}
//-----------------------------------------------------------------------------


short ISpinPattern::GetSpinB(int run, int Xing)
{
  return _SpinPat[run]->GetSpinB(Xing); 
}
//-----------------------------------------------------------------------------


short ISpinPattern::GetSpinY(int run, int Xing)
{
  return _SpinPat[run]->GetSpinY(Xing);
}
//-----------------------------------------------------------------------------


void ISpinPattern::GetSpinB(int run, short spin[XINGs])
{
  for (int i=0; i<XINGs; i++)
    spin[i] = _SpinPat[run]->GetSpinB(i);
}
//-----------------------------------------------------------------------------


void ISpinPattern::GetSpinY(int run, short spin[XINGs])
{
  for (int i=0; i<XINGs; i++)
    spin[i] = _SpinPat[run]->GetSpinY(i);
}
//-----------------------------------------------------------------------------


short ISpinPattern::GetNBunches(int run, SPIN B, SPIN Y)
{
  return _SpinPat[run]->GetNBunches(B, Y);
}
//-----------------------------------------------------------------------------


double ISpinPattern::GetBbcSum(SPIN B, SPIN Y)
{
  double bbcsum = 0.;

  MapSpinPat::const_iterator it;

  for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
    bbcsum += it->second->GetBbcSum(B, Y);

  return bbcsum;
}
//-----------------------------------------------------------------------------


double ISpinPattern::GetBbcSum(int run)
{
  return _SpinPat[run]->GetBbcSum();
}
//-----------------------------------------------------------------------------


double ISpinPattern::GetBbcSum(int run, SPIN B, SPIN Y)
{
  return _SpinPat[run]->GetBbcSum(B, Y);
}
//-----------------------------------------------------------------------------


double ISpinPattern::GetBbcSum(SPIN B1, SPIN Y1, SPIN B2, SPIN Y2)
{
  double bbcsum = 0.;

  MapSpinPat::const_iterator it;

  for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
    bbcsum += it->second->GetBbcSum(B1, Y1, B2, Y2);

  return bbcsum;
}
//-----------------------------------------------------------------------------


double ISpinPattern::GetBbcSum(int run, SPIN B1, SPIN Y1, SPIN B2, SPIN Y2)
{
  return _SpinPat[run]->GetBbcSum(B1, Y1, B2, Y2);
}
//-----------------------------------------------------------------------------


void ISpinPattern::GetRelLumi(RLUMI R_config, double &R, double &dR)
{
  double L1 = 0.;
  double L2 = 1.e10;

  switch (R_config)
    {
    case R_ALL:           // L(++,--)/L(+-,-+)
      L1 = GetBbcSum(BPOS, YPOS, BNEG, YNEG);
      L2 = GetBbcSum(BPOS, YNEG, BNEG, YPOS);
      break;

    case R_PVIOLATING1 :  // L(++)/L(--)
      L1 = GetBbcSum(BPOS, YPOS);
      L2 = GetBbcSum(BNEG, YNEG);
      break;

    case R_PVIOLATING2 :  // L(+-)/L(-+)
      L1 = GetBbcSum(BPOS, YNEG);
      L2 = GetBbcSum(BNEG, YPOS);
      break;

    case R_ANBLUE :       // L(Bu)/L(Bd)
      L1 = GetBbcSum(BUP, YUP, BUP, YDN);
      L2 = GetBbcSum(BDN, YUP, BDN, YDN);
      break;

    case R_ANYELLOW :     // L(Yu)/L(Yd)
      L1 = GetBbcSum(BUP, YUP, BDN, YUP);
      L2 = GetBbcSum(BUP, YDN, BDN, YDN);
      break;

    case R1_ANBLUE :      // L(Bu,Yu)/L(Bu,Yd)
    case R2_ANYELLOW :    // L(Bu,Yu)/L(Bd,Yu)
      L1 = GetBbcSum(BUP, YUP);
      L2 = GetBbcSum(BUP, YDN);
      break;

    case R2_ANBLUE :      // L(Bu,Yu)/L(Bd,Yu)
    case R1_ANYELLOW :    // L(Yu,Bu)/L(Yu,Bd)
      L1 = GetBbcSum(BUP, YUP);
      L2 = GetBbcSum(BDN, YUP);
      break;

    case R3_ANBLUE :      // L(Yu,Bu)/L(Yd,Bu)
    case R3_ANYELLOW :    // L(Yu,Bu)/L(Yd,Bd)
      L1 = GetBbcSum(BUP, YUP);
      L2 = GetBbcSum(BDN, YDN);
      break;

    default:
      L1 = 0.;
      L2 = 1.e10;
    }

  R = L1/L2;
  dR = R*sqrt(1/L1 + 1/L2);
}
//-----------------------------------------------------------------------------


void ISpinPattern::GetRelLumi(int run, RLUMI R_config, double &R, double &dR)
{
  double L1 = 0.;
  double L2 = 1.e10;

  switch (R_config)
    {
    case R_ALL:           // L(++,--)/L(+-,-+)
      L1 = _SpinPat[run]->GetBbcSum(BPOS, YPOS, BNEG, YNEG);
      L2 = _SpinPat[run]->GetBbcSum(BPOS, YNEG, BNEG, YPOS);
      break;

    case R_PVIOLATING1 :  // L(++)/L(--)
      L1 = _SpinPat[run]->GetBbcSum(BPOS, YPOS);
      L2 = _SpinPat[run]->GetBbcSum(BNEG, YNEG);
      break;

    case R_PVIOLATING2 :  // L(+-)/L(-+)
      L1 = _SpinPat[run]->GetBbcSum(BPOS, YNEG);
      L2 = _SpinPat[run]->GetBbcSum(BNEG, YPOS);
      break;

    case R_ANBLUE :       // L(Bu)/L(Bd)
      L1 = _SpinPat[run]->GetBbcSum(BUP, YUP, BUP, YDN);
      L2 = _SpinPat[run]->GetBbcSum(BDN, YUP, BDN, YDN);
      break;

    case R_ANYELLOW :     // L(Yu)/L(Yd)
      L1 = _SpinPat[run]->GetBbcSum(BUP, YUP, BDN, YUP);
      L2 = _SpinPat[run]->GetBbcSum(BUP, YDN, BDN, YDN);
      break;

    case R1_ANBLUE :      // L(Bu,Yu)/L(Bu,Yd)
    case R2_ANYELLOW :    // L(Bu,Yu)/L(Bd,Yu)
      L1 = _SpinPat[run]->GetBbcSum(BUP, YUP);
      L2 = _SpinPat[run]->GetBbcSum(BUP, YDN);
      break;

    case R2_ANBLUE :      // L(Bu,Yu)/L(Bd,Yu)
    case R1_ANYELLOW :    // L(Yu,Bu)/L(Yu,Bd)
      L1 = _SpinPat[run]->GetBbcSum(BUP, YUP);
      L2 = _SpinPat[run]->GetBbcSum(BDN, YUP);
      break;

    case R3_ANBLUE :      // L(Yu,Bu)/L(Yd,Bu)
    case R3_ANYELLOW :    // L(Yu,Bu)/L(Yd,Bd)
      L1 = _SpinPat[run]->GetBbcSum(BUP, YUP);
      L2 = _SpinPat[run]->GetBbcSum(BDN, YDN);
      break;

    default:
      L1 = 0.;
      L2 = 1.e10;
    }

  R = L1/L2;
  dR = R*sqrt(1/L1 + 1/L2);
}
//-----------------------------------------------------------------------------


void ISpinPattern::GetBbcIN(int run, float bbc[XINGs])
{
  for (int i=0; i<XINGs; i++)
    bbc[i] = _SpinPat[run]->GetBbcIn(i);
}
//-----------------------------------------------------------------------------


void ISpinPattern::GetPolarization(int run, float &PolB,
				   float &ePolBstat, float &ePolBsyst,
				   float &PolY,
				   float &ePolYstat, float &ePolYsyst)
{
  _SpinPat[run]->GetPolarization(PolB, ePolBstat, ePolBsyst,
				 PolY, ePolYstat, ePolYsyst);
}
//-----------------------------------------------------------------------------


void ISpinPattern::GetPolarizationB(int run, float &PolB,
				    float &ePolBstat, float &ePolBsyst)
{
  _SpinPat[run]->GetPolarizationB(PolB, ePolBstat, ePolBsyst);
}
//-----------------------------------------------------------------------------


void ISpinPattern::GetPolarizationY(int run, float &PolY,
				    float &ePolYstat, float &ePolYsyst)
{
  _SpinPat[run]->GetPolarizationY(PolY, ePolYstat, ePolYsyst);
}
//-----------------------------------------------------------------------------


void ISpinPattern::DrawBBC(int run)
{
  if (!Exists(run)) return;

  float xing[XINGs];
  float bbc[XINGs];

  for (int i=0; i<XINGs; i++)
    {
      bbc[i] = _SpinPat[run]->GetBbcIn(i);
      xing[i] = i*1.0;
    }

  TGraph *gr = new TGraph(XINGs, xing, bbc);
  TString title = "bbc_in, run/fill ";
  title += run;
  if (_Rndm) title += ", random";
  if (_Even) title += ", even bunches";
  if (_Odd ) title += ", odd bunches";

  if (gPad) gPad->Clear();

  gr->SetTitle(title.Data());
  gr->GetXaxis()->SetLimits(-1.,120.);
  gr->SetFillColor(43);
  gr->Draw("AB");
  gr->Draw("*");
}
//-----------------------------------------------------------------------------


void ISpinPattern::DrawSpinB(int run)
{
  if (!Exists(run)) return;

  int xing[XINGs];
  int spinB[XINGs];

  for (int i=0; i<XINGs; i++)
    {
      spinB[i] = static_cast<int>(_SpinPat[run]->GetSpinB(i));
      xing[i] = i;
    }

  TGraph *gr = new TGraph(XINGs, xing, spinB);
  TString title = "spin pattern blue, run ";
  title += run;
  if (_Rndm) title += ", random";
  if (_Even) title += ", even bunches";
  if (_Odd ) title += ", odd bunches";

  if (gPad) gPad->Clear();

  gr->SetTitle(title.Data());
  gr->GetXaxis()->SetLimits(-1.,120.);
  gr->SetFillColor(43);
  gr->Draw("AB");
  gr->Draw("*");
}
//-----------------------------------------------------------------------------


void ISpinPattern::DrawSpinY(int run)
{
  if (!Exists(run)) return;

  int xing[XINGs];
  int spinY[XINGs];

  for (int i=0; i<XINGs; i++)
    {
      spinY[i] = static_cast<int>(_SpinPat[run]->GetSpinY(i));
      xing[i] = i;
    }

  TGraph *gr = new TGraph(XINGs, xing, spinY);
  TString title = "spin pattern yellow, run ";
  title += run;
  if (_Rndm) title += ", random";
  if (_Even) title += ", even bunches";
  if (_Odd ) title += ", odd bunches";

  if (gPad) gPad->Clear();

  gr->SetTitle(title.Data());
  gr->GetXaxis()->SetLimits(-1.,120.);
  gr->SetFillColor(43);
  gr->Draw("AB");
  gr->Draw("*");
}
//-----------------------------------------------------------------------------


void ISpinPattern::Print(int run)
{
  if (!Exists(run)) return;
  _SpinPat[run]->Print();
}
//-----------------------------------------------------------------------------


TGraphErrors* ISpinPattern::RelLumiGraph(RLUMI R_config)
{
  int N = _SpinPat.size();

  double *relLumi  = new double[N];
  double *relLumiE = new double[N];
  double *run_fill = new double[N];

  MapSpinPat::const_iterator it;

  int index = 0;
  for (it = _SpinPat.begin(); it != _SpinPat.end(); it++)
    {
      run_fill[index] = float(it->first);
      GetRelLumi(it->first, R_config, relLumi[index], relLumiE[index]);

      index++;
    }

  TGraphErrors *gr = new TGraphErrors(index, run_fill, relLumi, 0, relLumiE);
  gr->SetMaximum(1.5);
  gr->SetMinimum(0.5);
  gr->SetLineColor(41);
  gr->SetMarkerColor(4);

  return gr;
}
//-----------------------------------------------------------------------------


void ISpinPattern::ReadGoodRunsList(const char* file)
{
  cout << "ISpinPattern::ReadGoodRunsList" << endl;
  _nGoodRuns = 0;
  ifstream fin(file);

  while (fin >> _GoodRun[_nGoodRuns])
    _nGoodRuns++;

  fin.close();
  cout << _nGoodRuns << " good runs!" << endl;
}
//-----------------------------------------------------------------------------


bool ISpinPattern::GoodRun(int run)
{
  bool good = false;

  for (int i=0; i<_nGoodRuns; i++)
    if (_GoodRun[i] == run)
      {
	good = true;
	break;
      }

  return good;
}
//-----------------------------------------------------------------------------
//=============================================================================
//-----------------------------------------------------------------------------


ISpinPatternRun::ISpinPatternRun(int run, short *spinB, short *spinY,
				 float *bbc_in, float *zdc_in) :
  _RunNum(run),
  _PolB(0),
  _PolY(0),
  _ePolBstat(0),
  _ePolBsyst(0),
  _ePolYstat(0),
  _ePolYsyst(0)
{
  _SpinBdb = new short[XINGs];
  _SpinYdb = new short[XINGs];
  _SpinB   = new short[XINGs];
  _SpinY   = new short[XINGs];
  _BbcIn   = new float[XINGs];
  _ZdcIn   = new float[XINGs];

  SetFillPattern(spinB, spinY);

  for (int i=0; i<XINGs; i++)
    {
      //if (spinB[i] != 1 && spinB[i] != -1) spinB[i] = 0;
      //if (spinY[i] != 1 && spinY[i] != -1) spinY[i] = 0;

      _SpinBdb[i] = spinB[i];
      _SpinYdb[i] = spinY[i];

      // if (bbc_in[i] > 0)
      // 	{
      // 	  _SpinBdb[i] = spinB[i];
      // 	  _SpinYdb[i] = spinY[i];
      // 	}
      // else
      // 	{
      // 	  _SpinBdb[i] = 0;
      // 	  _SpinYdb[i] = 0;
      // 	}

      _SpinB[i] = _SpinBdb[i];
      _SpinY[i] = _SpinYdb[i];

      if (_SpinBdb[i] == 0 || _SpinYdb[i] == 0)
	{
	  _BbcIn[i] = 0.;
	  _ZdcIn[i] = 0.;
	}
      else
	{
	  _BbcIn[i] = bbc_in[i];
	  _ZdcIn[i] = zdc_in[i];
	}
    }
}
//-----------------------------------------------------------------------------


ISpinPatternRun::~ISpinPatternRun()
{
  delete [] _SpinBdb;
  delete [] _SpinYdb;
  delete [] _SpinB;
  delete [] _SpinY;
  delete [] _BbcIn;
  delete [] _ZdcIn;
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::Reset()
{
  for (int i=0; i<XINGs; i++)
    {
      _SpinB[i] = _SpinBdb[i];
      _SpinY[i] = _SpinYdb[i];
    }
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::Shuffle(TRandom3 *rnd)
{
  for (int i=0; i<XINGs; i++)
    _SpinB[i] = (rnd->Rndm() > 0.5) ? 1 : -1;

  for (int i=0; i<XINGs; i++)
    _SpinY[i] = (rnd->Rndm() > 0.5) ? 1 : -1;
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::KnuthShuffle(TRandom3 *rnd)
{
  int bunch1[120], bunch2[120];
  int nFilledBunch = 0;

  for (int i=0; i<XINGs; i++)
    {
      _SpinB[i] = _SpinBdb[i];
      _SpinY[i] = _SpinYdb[i];

      if (_SpinBdb[i] != 0 && _SpinYdb[i] != 0)
	{
	  bunch1[nFilledBunch] = i;
	  bunch2[nFilledBunch] = i;
	  nFilledBunch++;
	}
    }

  for (int i=nFilledBunch-1; i>=0; i--)
    {
      int random = rnd->Integer(i);
      int temp = bunch1[i];
      bunch1[i] = bunch1[random];
      bunch1[random] = temp;
    }

  for (int i=0; i<nFilledBunch; i++)
    {
      _SpinB[bunch2[i]] = _SpinBdb[bunch1[i]];
      _SpinY[bunch2[i]] = _SpinYdb[bunch1[i]];
    }
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::SetEvenBunches()
{
  for (int i=0; i<XINGs; i++)
    if (i%2 == 0)
      {
	_SpinB[i] = _SpinBdb[i];
	_SpinY[i] = _SpinYdb[i];
      }
    else
      {
	_SpinB[i] = 0;
	_SpinY[i] = 0;
      }
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::SetOddBunches()
{
  for (int i=0; i<XINGs; i++)
    if (i%2 == 1)
      {
	_SpinB[i] = _SpinBdb[i];
	_SpinY[i] = _SpinYdb[i];
      }
    else
      {
	_SpinB[i] = 0;
	_SpinY[i] = 0;
      }
}
//-----------------------------------------------------------------------------


double ISpinPatternRun::GetBbcSum()
{
  double bbcsum = 0.;

  for (int i=0; i<XINGs; i++)
    bbcsum += double(_BbcIn[i]);

  return bbcsum;
}
//-----------------------------------------------------------------------------


double ISpinPatternRun::GetBbcSum(short B, short Y)
{
  double bbcsum = 0.;

  for (int i=0; i<XINGs; i++)
    if (_SpinB[i] == B && _SpinY[i] == Y)
      bbcsum += double(_BbcIn[i]);

  return bbcsum;
}
//-----------------------------------------------------------------------------


double ISpinPatternRun::GetBbcSum(short B1, short Y1, short B2, short Y2)
{
  return GetBbcSum(B1, Y1) + GetBbcSum(B2, Y2);
}
//-----------------------------------------------------------------------------


short ISpinPatternRun::GetNBunches(short B, short Y)
{
  short n = 0;
  for (int i=0; i<XINGs; i++)
    if(_SpinB[i] == B && _SpinY[i] == Y) n++;

  return n;
}
//-----------------------------------------------------------------------------


float ISpinPatternRun::GetBbcIn(int Xing)
{
  if (_SpinB[Xing] !=0 || _SpinY[Xing] !=0) return _BbcIn[Xing];
  else return 0.;
}
//-----------------------------------------------------------------------------


float ISpinPatternRun::GetZdcIn(int Xing)
{
  if (_SpinB[Xing] !=0 || _SpinY[Xing] !=0) return _ZdcIn[Xing];
  else return 0.;
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::SetFillPattern(short *spinB, short *spinY)
{
  _FillPat = -9;

  TString patB1 = "-+-++-+-";
  TString patB2 = "+-+--+-+";
  TString patB3 = "-+-+-+-+";
  TString patB4 = "+-+-+-+-";
  TString patY1 = "--++--++";
  TString patY2 = "++--++--";

  TString patB = "";
  TString patY = "";

  short k=0;
  for (int i=0; i<XINGs; i++)
    {
      if (spinB[i] == 0 || spinY[i] == 0) continue;

      if (k<8)
	{
	  if (spinB[i] == +1) patB += "+";
	  if (spinB[i] == -1) patB += "-";
	  if (spinY[i] == +1) patY += "+";
	  if (spinY[i] == -1) patY += "-";
	  k++;
	}
      else
	break;
    }

  if ( (patB.Contains(patB1) || patB.Contains(patB2) || // run6, run8
	patB.Contains(patB3) || patB.Contains(patB4)) && // run5
       (patY.Contains(patY1) || patY.Contains(patY2)) )
    {
      patB.Remove(4,4);
      patY.Remove(4,4);

      if (patB.Contains("-+-+") &&
	  patY.Contains("--++") ) _FillPat = 1;

      if (patB.Contains("+-+-") &&
	  patY.Contains("--++") ) _FillPat = 2;

      if (patB.Contains("-+-+") &&
	  patY.Contains("++--") ) _FillPat = 3;

      if (patB.Contains("+-+-") &&
	  patY.Contains("++--") ) _FillPat = 4;
    }
  else
    {
      patB = "";
      patY = "";

      for (int i=0; i<8; i++)
	{
	  if (spinB[i] == +1) patB += "+";
	  if (spinB[i] == -1) patB += "-";
	  if (spinB[i] ==  0) patB += ".";
	  if (spinY[i] == +1) patY += "+";
	  if (spinY[i] == -1) patY += "-";
	  if (spinY[i] ==  0) patY += ".";

	  if (spinB[i] == 0 || spinY[i] == 0)
	    {
	      patB1.Replace(i, 1, ".");
	      patB2.Replace(i, 1, ".");
	      patB3.Replace(i, 1, ".");
	      patB4.Replace(i, 1, ".");
	      patY1.Replace(i, 1, ".");
	      patY2.Replace(i, 1, ".");
	    }
	}

      if ( (patB.Contains(patB1) && patY.Contains(patY1)) ||
	   (patB.Contains(patB3) && patY.Contains(patY1)) )
	_FillPat = 1;

      if ( (patB.Contains(patB2) && patY.Contains(patY1)) ||
	   (patB.Contains(patB4) && patY.Contains(patY1)) )
	_FillPat = 2;

      if ( (patB.Contains(patB1) && patY.Contains(patY2)) ||
	   (patB.Contains(patB3) && patY.Contains(patY2)) )
	_FillPat = 3;

      if ( (patB.Contains(patB2) && patY.Contains(patY2)) ||
	   (patB.Contains(patB4) && patY.Contains(patY2)) )
	_FillPat = 4;
    }
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::Print()
{
  TString spinB = "B  ";
  TString spinY = "Y  ";

  int same=0, diff=0, bup=0, bdn=0, yup=0, ydn=0;

  for (int i=0; i<XINGs; i++)
    {
      int polB, polY;
      polB = _SpinB[i];
      polY = _SpinY[i];

      if (polB*polY == +1) same++;
      if (polB*polY == -1) diff++;

      if (polB == +1) bup++;
      if (polB == -1) bdn++;
      if (polY == +1) yup++;
      if (polY == -1) ydn++;

      if (polB == +1) spinB += "+";
      if (polB == -1) spinB += "-";
      if (polB ==  0) spinB += ".";

      if (polY == +1) spinY += "+";
      if (polY == -1) spinY += "-";
      if (polY ==  0) spinY += ".";
    }

  //cout << endl;
  cout << "spin pattern for run/fill " << _RunNum 
       << "; Fill pattern " << _FillPat << endl;
  cout << spinB.Data() << endl;
  cout << spinY.Data() << endl;
  cout << endl;

  cout << "same helicity bunches = " << same << endl;
  cout << "diff helicity bunches = " << diff << endl;
  cout << "No. of B+ bunches = " << bup << endl;
  cout << "No. of B- bunches = " << bdn << endl;
  cout << "No. of Y+ bunches = " << yup << endl;
  cout << "No. of Y- bunches = " << ydn << endl;
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::SetPolarization(float polB,
				      float ePolBstat, float ePolBsyst,
				      float polY,
				      float ePolYstat, float ePolYsyst)
{
  _PolB      = polB;
  _ePolBstat = ePolBstat;
  _ePolBsyst = ePolBsyst;
  _PolY      = polY;
  _ePolYstat = ePolYstat;
  _ePolYsyst = ePolYsyst;
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::GetPolarization(float &PolB,
				      float &ePolBstat, float &ePolBsyst,
				      float &PolY,
				      float &ePolYstat, float &ePolYsyst)
{
  PolB      = _PolB;
  ePolBstat = _ePolBstat;
  ePolBsyst = _ePolBsyst;
  PolY      = _PolY;
  ePolYstat = _ePolYstat;
  ePolYsyst = _ePolYsyst;
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::GetPolarizationB(float &PolB,
				       float &ePolBstat, float &ePolBsyst)
{
  PolB      = _PolB;
  ePolBstat = _ePolBstat;
  ePolBsyst = _ePolBsyst;
}
//-----------------------------------------------------------------------------


void ISpinPatternRun::GetPolarizationY(float &PolY,
				       float &ePolYstat, float &ePolYsyst)
{
  PolY      = _PolY;
  ePolYstat = _ePolYstat;
  ePolYsyst = _ePolYsyst;
}
//-----------------------------------------------------------------------------
