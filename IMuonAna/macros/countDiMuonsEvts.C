#include <iostream>
#include <fstream>
#include <math.h>


#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

using namespace std;

#define MINMS  2.2
#define MAXMS  4.8

enum PAR {A, K, Ajpsi, Mjpsi, Sjpsi, Apsi, Mpsi, Spsi};
enum PAR2 {P0, P1, P2, P3};


double GausExpoFunc(double *x, double *p);
double GausPol3Func(double *x, double *p);


void countDiMuonsEvts()
{
  int    RunNumber;
  int    FillNumber;
  int    EvtNumber;
  short  NMuons ;
  short  SpinXingID;
  short  SpinB;
  short  SpinY;
  double Zvtx;

  double Mu1_DDG0;
  double Mu1_DG0;
  double Mu1_DS3;
  double Mu1_DS3ctp;
  double Mu1_MuTrChi2;
  double Mu1_MuIdChi2;
  double Mu1_Px;
  double Mu1_Py;
  double Mu1_Pz;
  double Mu1_Pt;
  double Mu1_P;
  short  Mu1_Charge;
  short  Mu1_nMuTrHits;
  short  Mu1_nMuIdHits;
  unsigned int Mu1_MuTrHitPat;
  unsigned int Mu1_MuIdHitPat;

  double Mu2_DDG0;
  double Mu2_DG0;
  double Mu2_DS3;
  double Mu2_DS3ctp;
  double Mu2_MuTrChi2;
  double Mu2_MuIdChi2;
  double Mu2_Px;
  double Mu2_Py;
  double Mu2_Pz;
  double Mu2_Pt;
  double Mu2_P;
  short  Mu2_Charge;
  short  Mu2_nMuTrHits;
  short  Mu2_nMuIdHits;
  unsigned int Mu2_MuTrHitPat;
  unsigned int Mu2_MuIdHitPat;

  short  diMuCharge;
  double diMuMass;
  double diMuP;
  double diMuPt;
  double diMuPz;
  double diMuPhi;
  double diMuEta;


  ifstream fin("junk.txt");
  float *runNum = new float[300];

  int runCount = 0;
  while (fin >> runNum[runCount]) runCount++;

  float *nDiMu = new float[runCount];

  for (int i=0; i<runCount; i++)
    nDiMu[i] = 0.;

  TH1D *hMass = new TH1D("hMass", "hMass",100,1.,6);


  // TF1 *fMass  = new TF1("fMassUNM",  GausExpoFunc, MINMS, MAXMS, 6);

  // fMass->SetParNames("A","K","Njpsi","Mjpsi","Sjpsi",
  // 			"Apsi","Mpsi","Spsi");

  // fMass->SetParameter(A, 100.);
  // fMass->SetParameter(K, 1.);

  // fMass->SetParameter(Ajpsi, 200.);
  // fMass->SetParameter(Mjpsi, 3.1);
  // fMass->SetParameter(Sjpsi, 0.12);

  // fMass->SetParameter(Apsi, 1.);
  // fMass->SetParLimits(Apsi, 0., 9999999.);


  TF1 *fMass  = new TF1("fMass", GausPol3Func, MINMS, MAXMS, 8);

  fMass->SetParNames("P0","P1","P2","P3","Njpsi","Mjpsi","Sjpsi",
			"Npsi","Mpsi","Spsi");

  fMass->SetParameter(P0, 100.);
  fMass->SetParameter(P1, -10.);
  fMass->SetParameter(P2, 1.);
  fMass->SetParameter(P3, 1.);

  fMass->SetParameter(Ajpsi+2, 6000.);
  fMass->SetParameter(Mjpsi+2, 3.1);
  fMass->SetParameter(Sjpsi+2, 0.12);

  fMass->SetParameter(Apsi+2, 60.);
  fMass->SetParLimits(Apsi+2, 0., 9999999.);


  //TFile *f = new TFile(fileName[nfile]);
  //TTree *jpsi = (TTree*)f->Get("dimuons");

  TChain *jpsi = new TChain("dimuons");
  jpsi->Add("diMuonAN_taxi105_1.root");
  jpsi->Add("diMuonAN_taxi105_2.root");


  jpsi->SetBranchAddress("RunNumber",  &RunNumber);
  jpsi->SetBranchAddress("FillNumber", &FillNumber);
  jpsi->SetBranchAddress("EvtNumber",  &EvtNumber);
  jpsi->SetBranchAddress("NMuons",     &NMuons);
  jpsi->SetBranchAddress("SpinXingIS", &SpinXingID);
  jpsi->SetBranchAddress("SpinB",      &SpinB);
  jpsi->SetBranchAddress("SpinY",      &SpinY);
  jpsi->SetBranchAddress("Zvtx",       &Zvtx);

  jpsi->SetBranchAddress("Mu1_DDG0",     &Mu1_DDG0);
  jpsi->SetBranchAddress("Mu1_DG0",      &Mu1_DG0);
  jpsi->SetBranchAddress("Mu1_DS3",      &Mu1_DS3);
  jpsi->SetBranchAddress("Mu1_DS3ctp",   &Mu1_DS3ctp);
  jpsi->SetBranchAddress("Mu1_MuTrChi2", &Mu1_MuTrChi2);
  jpsi->SetBranchAddress("Mu1_MuIdChi2", &Mu1_MuIdChi2);
  jpsi->SetBranchAddress("Mu1_Px",       &Mu1_Px);
  jpsi->SetBranchAddress("Mu1_Py",       &Mu1_Py);
  jpsi->SetBranchAddress("Mu1_Pz",       &Mu1_Pz);
  jpsi->SetBranchAddress("Mu1_Pt",       &Mu1_Pt);
  jpsi->SetBranchAddress("Mu1_P",        &Mu1_P);
  jpsi->SetBranchAddress("Mu1_Charge",   &Mu1_Charge);

  jpsi->SetBranchAddress("Mu1_nMuTrHits",  &Mu1_nMuTrHits);
  jpsi->SetBranchAddress("Mu1_nMuIdHits",  &Mu1_nMuIdHits);
  jpsi->SetBranchAddress("Mu1_MuTrHitPat", &Mu1_MuTrHitPat);
  jpsi->SetBranchAddress("Mu1_MuIdHitPat", &Mu1_MuIdHitPat);

  jpsi->SetBranchAddress("Mu2_DDG0",     &Mu2_DDG0);
  jpsi->SetBranchAddress("Mu2_DG0",      &Mu2_DG0);
  jpsi->SetBranchAddress("Mu2_DS3",      &Mu2_DS3);
  jpsi->SetBranchAddress("Mu2_DS3ctp",   &Mu2_DS3ctp);
  jpsi->SetBranchAddress("Mu2_MuTrChi2", &Mu2_MuTrChi2);
  jpsi->SetBranchAddress("Mu2_MuIdChi2", &Mu2_MuIdChi2);
  jpsi->SetBranchAddress("Mu2_Px",       &Mu2_Px);
  jpsi->SetBranchAddress("Mu2_Py",       &Mu2_Py);
  jpsi->SetBranchAddress("Mu2_Pz",       &Mu2_Pz);
  jpsi->SetBranchAddress("Mu2_Pt",       &Mu2_Pt);
  jpsi->SetBranchAddress("Mu2_P",        &Mu2_P);
  jpsi->SetBranchAddress("Mu2_Charge",   &Mu2_Charge);

  jpsi->SetBranchAddress("Mu2_nMuTrHits",  &Mu2_nMuTrHits);
  jpsi->SetBranchAddress("Mu2_nMuIdHits",  &Mu2_nMuIdHits);
  jpsi->SetBranchAddress("Mu2_MuTrHitPat", &Mu2_MuTrHitPat);
  jpsi->SetBranchAddress("Mu2_MuIdHitPat", &Mu2_MuIdHitPat);

  jpsi->SetBranchAddress("diMuCharge",  &diMuCharge);
  jpsi->SetBranchAddress("diMuMass",    &diMuMass);
  jpsi->SetBranchAddress("diMuP",       &diMuP);
  jpsi->SetBranchAddress("diMuPt",      &diMuPt);
  jpsi->SetBranchAddress("diMuPz",      &diMuPz);
  jpsi->SetBranchAddress("diMuEta",     &diMuEta);

  cout << jpsi->GetEntries() << endl;


  int nRuns=0;
  int oldRun=-9;
  bool runExists = false;

  for (int evt=0; evt<jpsi->GetEntries(); evt++)
    {
      jpsi->GetEntry(evt);

      if (evt == 0)
	{
	  oldRun = RunNumber;
	}

      runExists = false;
      for (int i=0; i<runCount; i++)
	if (int(runNum[i]) == RunNumber) runExists = true;

      if (!runExists) continue;

      if (fabs(Zvtx) > 35.) continue;
      if (diMuCharge != 0) continue;
      if (Mu1_Pz*Mu2_Pz<0) continue;

      if (fabs(diMuEta)<1.2 || fabs(diMuEta)>2.2) continue;
      if (fabs(Zvtx) == 0 ) continue;

      if (diMuPt == 0  ||  diMuPt > 6) continue;
      if (diMuP > 80 )              continue;
      if (fabs(Mu1_Pz) < 1.4  || fabs(Mu2_Pz) < 1.4)    continue;
      if (fabs(Mu1_Pz) > 20   || fabs(Mu2_Pz) > 20 )    continue;
      if (Mu1_MuIdHitPat < 64 || Mu2_MuIdHitPat < 64 )  continue;
      if (Mu1_DDG0 > 10       || Mu2_DDG0 > 10)         continue;
      if (Mu1_MuTrChi2 > 30   || Mu2_MuTrChi2 > 30)     continue;

      if (Mu1_Pz < 0) continue;


      if (Mu1_Pz>0 && (fabs(Mu1_DG0) > 25 || fabs(Mu2_DG0) > 25))
      	continue;

      if (Mu1_Pz<0 && (fabs(Mu1_DG0) > 30 || fabs(Mu2_DG0) > 30))
      	continue;

      if (diMuMass < 1.5 || diMuMass > 6.) continue;

      nDiMu[nRuns]++;

      if (diMuMass > 1. && diMuMass < 6.) hMass->Fill(diMuMass);

      if (oldRun != RunNumber)
      	{
      	  //cout << oldRun << endl;
      	  if (nDiMu[nRuns] == 0) cout << oldRun << endl;
      	  oldRun = RunNumber;
      	  nRuns++;
      	}

    }
  nRuns++;

  cout << nRuns << "   " << runCount << endl;


  TGraph *gN = new TGraph(nRuns, runNum, nDiMu);

  TCanvas *cc = new TCanvas("cc","cc",1300,700);
  cc->Divide(2,2);

  cc->cd(1);
  gN->Draw("AP");
  cc->cd(2);
  hMass->Draw();
  hMass->Fit(fMass,"R");

}
//-----------------------------------------------------------------------------


double GausExpoFunc(double *m, double *par)
{
  double x = m[0];

  double A = par[0];
  double K = par[1];

  double N1 = par[2];
  double M1 = par[3];
  double S1 = par[4];

  double N2 = par[5];
  double M2 = M1*(3.686/3.097);
  double S2 = (M2/M1)*S1;

  double f = A*exp(-K*x) +
    N1*(5./100.)*TMath::Gaus(x, M1, S1, true) +
    N2*(5./100.)*TMath::Gaus(x, M2, S2, true) ;

  return f;
}
//-----------------------------------------------------------------------------


double GausPol3Func(double *m, double *par)
{
  double x = m[0];

  double P0 = par[0];
  double P1 = par[1];
  double P2 = par[2];
  double P3 = par[3];

  double N1 = par[4];
  double M1 = par[5];
  double S1 = par[6];

  double N2 = par[7];
  double M2 = M1*(3.686/3.097);
  double S2 = (M2/M1)*S1;

  double f = P0 + P1*x + P2*x*x + P3*x*x*x +
    0.05*N1*TMath::Gaus(x, M1, S1, true) +
    0.05*N2*TMath::Gaus(x, M2, S2, true) ;

  return f;
}
//-----------------------------------------------------------------------------
