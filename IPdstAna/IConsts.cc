#include <iostream>
#include <stdio.h>
#include <string.h>


#include "IConsts.hh"
#include "TString.h"


using namespace std;

int     IConsts::RUN        =  0                      ;
int     IConsts::NPtTBins   =  8                      ;
int     IConsts::NPtABins   =  1                      ;
double* IConsts::PtT        =  0                      ;
double* IConsts::PtA        =  0                      ;
double  IConsts::MinPhi     = -PIby2                  ;
double  IConsts::MaxPhi     =  3.0*PI/2.0             ;
double  IConsts::MinDPhi    = -PIby2                  ;
double  IConsts::MaxDPhi    =  3.0*PI/2.0             ;
double  IConsts::MinDPhiN   = -PIby2                  ;
double  IConsts::MaxDPhiN   =  PIby2                  ;
double  IConsts::MinDPhiF   =  PIby2                  ;
double  IConsts::MaxDPhiF   =  3.0*PI/2.0             ;
double  IConsts::MinTheta   =  65.0                   ;
double  IConsts::MaxTheta   =  115.0                  ;
double  IConsts::EMCB1      = -33.75*PI/180.0         ;
double  IConsts::EMCB2      =  56.25*PI/180.0         ;
double  IConsts::DCHB1      = -38.75*PI/180.0         ;
double  IConsts::DCHB2      =  61.25*PI/180.0         ;
double  IConsts::NSigN      =  3.0                    ;
double  IConsts::NSigF      =  3.0                    ;
double* IConsts::SigN       =  0                      ;
double* IConsts::SigF       =  0                      ;
double  IConsts::EAsym      =  0.8                    ;
double  IConsts::Pi0Width   =  2.0                    ;
double  IConsts::ConeAngle  =  0.5                    ;
double  IConsts::MatchSig   =  3.0                    ;
short   IConsts::TrChar     =  0                      ;
short   IConsts::EMCArm     =  0                      ;
double  IConsts::BbcChLo    =  0.0                    ;
double  IConsts::BbcChHi    =  1.e9                   ;
short   IConsts::BbcNLo     =  0                      ;
short   IConsts::BbcNHi     =  9999                   ;
short   IConsts::AsymConfig =  0                      ;
char*   IConsts::FilePrefix                           ;
char*   IConsts::RunList                              ;
char*   IConsts::InDir                                ;
char*   IConsts::OutDir                               ;
char*   IConsts::OutFile                              ;
char*   IConsts::SpinFile                             ;
char*   IConsts::PbScStats                            ;
char*   IConsts::PbGlStats                            ;
char*   IConsts::EMCWarnMap                           ;
char*   IConsts::FillTable                            ;


short   IConsts::H[10][4]                             ;
bool    IConsts::Defined    =  false                  ;

IConsts::IConsts(int run)
{
  RUN = run;

  cout << "IConsts::IConsts" << endl;
  cout << "Analysing Run " << RUN << endl;

  double ptt[] = {2.0, 2.5, 3.0, 3.5, 4.2, 5.2, 6.4, 8.0, 10.0};
  double pta[] = {0.2, 50.0};

  PtT = new double[NPtTBins+1];
  PtA = new double[NPtABins+1];

  for (int i=0; i<=NPtTBins; i++)
    PtT[i] = ptt[i];

  for (int i=0; i<=NPtABins; i++)
    PtA[i] = pta[i];


  SigN = new double[NPtTBins];
  SigF = new double[NPtTBins];


//   //string dummy = "goodRunsList.txt";
//   //RUNLIST    = dummy.c_str();

  FilePrefix = new char[100];
  RunList    = new char[100];
  InDir      = new char[100];
  OutDir     = new char[100];
  OutFile    = new char[100];
  SpinFile   = new char[100];
  PbScStats  = new char[100];
  PbGlStats  = new char[100];
  EMCWarnMap = new char[100];
  FillTable  = new char[100];

  strcpy(FilePrefix, "IPdstReco_");
  strcpy(RunList,    "goodRunsList.txt");
  strcpy(InDir,      "pdstRootFiles/");
  strcpy(OutDir,     "./");
  strcpy(OutFile,    "output.root");
  strcpy(SpinFile,   "GL1P_dataFile.root");
  strcpy(PbScStats,  "ERT_pi0Stats_PbSc.root");
  strcpy(PbGlStats,  "ERT_pi0Stats_PbGl.root");
  strcpy(EMCWarnMap, "EMCal_WarnMap.dat");
  strcpy(FillTable,  "/home/iyounus/workdir/analysis/unm_pdst_ana/IPdstAna/"
	 "IFillLookupTable.dat");

  for (int i=0; i<NPtTBins; i++)
    {
      //SigN[i] = sigN[i];
      //SigF[i] = sigF[i];
      SigN[i] = PIby2/NSigN;//sigN[i];
      SigF[i] = PIby2/NSigF;//sigF[i];
    }


  //                   B & Y | B & Y
  short h[10][4] = {{ +1, +1, -1, -1 },   // same helicity
		    { +1, -1, -1, +1 },   // opposite helicity
		    { +1, +1, +1, +1 },   // ++
		    { -1, -1, -1, -1 },   // --
		    { +1, -1, +1, -1 },   // +-
		    { -1, +1, -1, +1 },   // -+
		    { +1, +1, +1, -1 },   // B+
		    { -1, +1, -1, -1 },   // B-
		    { +1, +1, -1, +1 },   // Y+
		    { +1, -1, -1, -1 }};  // Y-


  for (int i=0; i<10; i++)
    for (int j=0; j<4; j++)
      IConsts::H[i][j] = h[i][j];

  Defined = true;
}
//-----------------------------------------------------------------------------


IConsts::~IConsts()
{
  if (PtT ) delete [] PtT;
  if (PtA ) delete [] PtA;
  if (SigN) delete [] SigN;
  if (SigF) delete [] SigF;
}
//-----------------------------------------------------------------------------


void IConsts::SetRunList(char* list)
{
  RunList = list;
}
//-----------------------------------------------------------------------------


void IConsts::SetInDir(char* dir)
{
  InDir = dir;
}
//-----------------------------------------------------------------------------


void IConsts::SetOutDir(char* dir)
{
  OutDir = dir;
}
//-----------------------------------------------------------------------------


void IConsts::SetOutFile(char* file)
{
  OutFile = file;
}
//-----------------------------------------------------------------------------


void IConsts::SetSpinFile(char *spin)
{
  SpinFile = spin;
}
//-----------------------------------------------------------------------------


void IConsts::SetPbScStats(char* stats)
{
  PbScStats = stats;
}
//-----------------------------------------------------------------------------


void IConsts::SetPbGlStats(char* stats)
{
  PbGlStats = stats;
}
//-----------------------------------------------------------------------------


void IConsts::SetEMCWarnMap(char* map)
{
  EMCWarnMap = map;
}
//-----------------------------------------------------------------------------


void IConsts::SetFillTable(char* fill)
{
  FillTable = fill;
}
//-----------------------------------------------------------------------------


void IConsts::SetFilePrefix(char* pre)
{
  FilePrefix = pre;
}
//-----------------------------------------------------------------------------


void IConsts::SetTriggerPtBins(int nbins, double *ptt)
{
  NPtTBins = nbins;
  if (PtT) delete [] PtT;
  PtT = new double[NPtTBins+1];
  for (int i=0; i<=NPtTBins; i++)
    PtT[i] = ptt[i];
}
//-----------------------------------------------------------------------------


void IConsts::SetAssociatedPtBins(int nbins, double *pta)
{
  NPtABins = nbins;
  if (PtA) delete [] PtA;
  PtA = new double[NPtABins+1];
  for (int i=0; i<=NPtABins; i++)
    PtA[i] = pta[i];
}
//-----------------------------------------------------------------------------


void IConsts::Print()
{
  cout << "\n**************************************************************" 
       << "\n- Analysing Run       " << RUN
       << "\n- File prefix         " << FilePrefix
       << "\n- Run list            " << RunList
       << "\n- Input dir           " << InDir
       << "\n- Output dir          " << OutDir
       << "\n- Spin info file      " << SpinFile
       << "\n- Pi0 info PbSc       " << PbScStats
       << "\n- Pi0 info PbGl       " << PbGlStats
       << "\n- No of PtT bins      " << NPtTBins
       << "\n- PtT                 ";

  for (int i=0; i<=NPtTBins; i++)
    cout << PtT[i] << ", ";

  cout << "\n- No of PtA bins      " << NPtABins
       << "\n- PtA                 ";

  for (int i=0; i<=NPtABins; i++)
    cout << PtA[i] << ", ";

  cout << "\n- Min Phi             " << MinPhi*180./PI
       << "\n- Max Phi             " << MaxPhi*180./PI
       << "\n- Min dPhi            " << MinDPhi*180./PI
       << "\n- Max dPhi            " << MaxDPhi*180./PI
       << "\n- Min thate           " << MinTheta
       << "\n- Max thate           " << MaxTheta
       << "\n- Track Charge        " << TrChar
       << "\n- EMC Arm             " << EMCArm
       << "\n- En asym             " << EAsym
       << "\n- Pi0 width           " << Pi0Width << " sigma"
       << "\n- Cone Angle          " << ConeAngle
       << "\n- PC3/EMC matching    " << MatchSig << " sigma"
       << "\n- Near side peak      " << NSigN    << " sigma"
       << "\n- Away side peak      " << NSigF    << " sigma"
       << "\n- Asym config         " << AsymConfig
       << "\n- EMC bound.1         " << EMCB1*180./PI
       << "\n- EMC bound.2         " << EMCB2*180./PI
       << "\n- DCH bound.1         " << DCHB1*180./PI
       << "\n- DCH bound.2         " << DCHB2*180./PI
       << "\n- Min Bbc char        " << BbcChLo
       << "\n- Max Bbc char        " << BbcChHi
       << "\n- Min Bbc mult        " << BbcNLo
       << "\n- Max Bbc mult        " << BbcNHi
       << "\n**************************************************************\n" 
       << endl;
}
//-----------------------------------------------------------------------------
