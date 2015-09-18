#include "IJpsiAN.hh"

#include <iostream>
#include <fstream>
#include <assert.h>


#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraphErrors.h"


#include "IConsts.hh"
#include "ISpinPattern.hh"


using namespace std;

double GausExpoFunc(double *x, double *p);


#define MINMS  2.2
#define MAXMS  4.8

#define BlueSpinShift  0.   // shift in spin direction from verticle
#define YeloSpinShift  0.   // shift in spin direction from verticle


IJpsiAN::IJpsiAN()
{
  assert(IConsts::Defined); // IConsts must be defined before calling this

  _NPtBins = IConsts::NPtTBins;

  _file = new TFile("diMuonAN_run8pp_JpsiMassRange.root");


  CreateArrays();
  GetHistos();
  SetMassFunc();
  CountJpsi();    // this should come after SetMassFunc().
  GetSpinInfo();  // this should come after CountJpsi(), because the fill list
                  // exists only after CountJpsi() is called.

  CalculateAccFunc();
  CalculateAsym();

}
//-----------------------------------------------------------------------------


void IJpsiAN::CreateArrays()
{
  cout << "IJpsiAN::CreateArrays" << endl;

  _Fill = new double[100]; // this is to plot graphs
  _nFills = 100;   // I don't expect the number fills to large than 100

  _Njpsi = new double**[_NPtBins];
  _Nbkgr = new double**[_NPtBins];

  _AN  = new double**[_NPtBins];
  _dAN = new double**[_NPtBins];

  _SumPhi    = new double**[_NPtBins];
  _SumCosPhi = new double**[_NPtBins];
  _AccFunc   = new double**[_NPtBins];
  _dAccFunc  = new double**[_NPtBins];


  for (int i=0; i<_NPtBins; i++)
    {
      _Njpsi[i] = new double*[nSPINCONF];
      _Nbkgr[i] = new double*[nSPINCONF];
      for (int j=0; j<nSPINCONF; j++)
	{
	  _Njpsi[i][j] = new double[_nFills];
	  _Nbkgr[i][j] = new double[_nFills];

	  for (int k=0; k<_nFills; k++)
	    {
	      _Njpsi[i][j][k] = 0.;
	      _Nbkgr[i][j][k] = 0.;
	    }
	}


      _AN[i]  = new double*[nASYMCONF2];
      _dAN[i] = new double*[nASYMCONF2];
      for (int j=0; j<nASYMCONF2; j++)
	{
	  _AN[i][j]  = new double[_nFills];
	  _dAN[i][j] = new double[_nFills];

	  for (int k=0; k<_nFills; k++)
	    {
	      _AN[i][j][k]  = 0.;
	      _dAN[i][j][k] = 0.;
	    }
	}

      _SumPhi[i]    = new double*[nASYMCONF1];
      _SumCosPhi[i] = new double*[nASYMCONF1];
      _AccFunc[i]   = new double*[nASYMCONF1];
      _dAccFunc[i]  = new double*[nASYMCONF1];
      for (int j=0; j<nASYMCONF1; j++)
	{
	  _SumPhi[i][j]    = new double[_nFills];
	  _SumCosPhi[i][j] = new double[_nFills];
	  _AccFunc[i][j]   = new double[_nFills];
	  _dAccFunc[i][j]  = new double[_nFills];

	  for (int k=0; k<_nFills; k++)
	    {
	      _SumPhi[i][j][k]    = 0.;
	      _SumCosPhi[i][j][k] = 0.;
	      _AccFunc[i][j][k]   = 0.;
	      _dAccFunc[i][j][k]  = 0.;
	    }
	}
    }
}
//-----------------------------------------------------------------------------


void IJpsiAN::GetSpinInfo()
{
  cout << "IJpsiAN::GetSpinInfo" << endl;
  _spin = new ISpinPattern(IConsts::SpinFile, IConsts::RunList, "Fill");

  for (int i=0; i<nBEAM; i++)
    {
      _Pol[i]      = new float [_nFills];
      _ePolStat[i] = new float [_nFills];
      _ePolSyst[i] = new float [_nFills];
      _R[i]        = new double[_nFills];
      _dR[i]       = new double[_nFills];
    }


  for (int i=0; i<_nFills; i++)
    {
      int fill = int(_Fill[i]);

      if (!_spin->Exists(fill))
	{
	  cout << fill << "  ************************" << endl;
	  continue;
	}

      _spin->
	GetPolarization(fill,
			_Pol[Blue][i], _ePolStat[Blue][i], _ePolSyst[Blue][i],
			_Pol[Yelo][i], _ePolStat[Yelo][i], _ePolSyst[Yelo][i]);

      _spin->
	GetRelLumi(fill, ISpinPattern::R_ANBLUE, _R[Blue][i], _dR[Blue][i]);

      _spin->
	GetRelLumi(fill, ISpinPattern::R_ANYELLOW, _R[Yelo][i], _dR[Yelo][i]);
    }
}
//-----------------------------------------------------------------------------


void IJpsiAN::GetHistos()
{
  cout << "IJpsiAN::GetHistos" << endl;

  _hJpsiMs[North] = new TH1D("hJpsiMs_N", "hJpsiMs_N", 100, 1., 6.);
  _hJpsiMs[South] = new TH1D("hJpsiMs_S", "hJpsiMs_S", 100, 1., 6.);

  _hDiMuMs[North] = new TH1D("hDiMuMs_N", "hDiMuMs_N", 100, 1., 6.);
  _hDiMuMs[South] = new TH1D("hDiMuMs_S", "hDiMuMs_S", 100, 1., 6.);

  // get all histograms from the file.
  TString hname = "";

  hname = "hJpsiMsN";
  _hJpsiMs[North]->Add((TH1D*)_file->Get(hname.Data()));

  hname = "hJpsiMsS";
  _hJpsiMs[South]->Add((TH1D*)_file->Get(hname.Data()));

  hname = "hDiMuMsN";
  _hDiMuMs[North]->Add((TH1D*)_file->Get(hname.Data()));

  hname = "hDiMuMsS";
  _hDiMuMs[South]->Add((TH1D*)_file->Get(hname.Data()));

}
//-----------------------------------------------------------------------------


void IJpsiAN::SetMassFunc()
{
  _fMs[North] = new TF1("fMsN", GausExpoFunc, MINMS, MAXMS, 6);
  _fMs[South] = new TF1("fMsS", GausExpoFunc, MINMS, MAXMS, 6);

  for (int i=0; i<nARM; i++)
    {
      _fMs[i]->SetParNames("A","K","Njpsi","Mjpsi","Sjpsi",
			   "Npsi","Mpsi","Spsi");

      _fMs[i]->SetParameter(A, 100.);
      _fMs[i]->SetParameter(K, 1.);

      _fMs[i]->SetParameter(Ajpsi, 200.);
      _fMs[i]->SetParameter(Mjpsi, 3.1);
      _fMs[i]->SetParameter(Sjpsi, 0.12);

      _fMs[i]->SetParameter(Apsi, 1.);
      _fMs[i]->SetParLimits(Apsi, 0., 9999999.);
    }


  TCanvas *cc = new TCanvas("cc","cc",1000,500);
  cc->Divide(2,1);

  for (int i=0; i<nARM; i++)
    {
      cc->cd(1+i);
      _hJpsiMs[i]->Draw();
      _hJpsiMs[i]->Fit(_fMs[i],"R");

      _JpsiMs[i]  = _fMs[i]->GetParameter(Mjpsi);
      _JpsiSig[i] = _fMs[i]->GetParameter(Sjpsi);
    }
}
//-----------------------------------------------------------------------------


void IJpsiAN::CountJpsi()
{
  cout << "IJpsiAN::CountJpsi" << endl;

  int    RunNumber;
  int    FillNumber;
  int    EvtNumber;
  short  NMuons ;
  short  SpinXingID;
  short  SpinB;
  short  SpinY;
  double Zvtx;

  short  diMuCharge;
  double diMuMass;
  double diMuP;
  double diMuPt;
  double diMuPz;
  double diMuPhi;
  double diMuEta;

  TTree *jpsi = (TTree *)_file->Get("dimuons");

  jpsi->SetBranchAddress("RunNumber",  &RunNumber);
  jpsi->SetBranchAddress("FillNumber", &FillNumber);
  jpsi->SetBranchAddress("EvtNumber",  &EvtNumber);
  jpsi->SetBranchAddress("NMuons",     &NMuons);
  jpsi->SetBranchAddress("SpinXingID", &SpinXingID);
  jpsi->SetBranchAddress("SpinB",      &SpinB);
  jpsi->SetBranchAddress("SpinY",      &SpinY);
  jpsi->SetBranchAddress("Zvtx",       &Zvtx);

  jpsi->SetBranchAddress("diMuCharge", &diMuCharge);
  jpsi->SetBranchAddress("diMuMass",   &diMuMass);
  jpsi->SetBranchAddress("diMuP",      &diMuP);
  jpsi->SetBranchAddress("diMuPt",     &diMuPt);
  jpsi->SetBranchAddress("diMuPz",     &diMuPz);
  jpsi->SetBranchAddress("diMuPhi",    &diMuPhi);
  jpsi->SetBranchAddress("diMuEta",    &diMuEta);


  double minMs[nARM];
  double maxMs[nARM];

  for (int i=0; i<nARM; i++)
    {
      minMs[i] = _JpsiMs[i] - 2.*_JpsiSig[i];
      maxMs[i] = _JpsiMs[i] + 2.*_JpsiSig[i];
    }


  short Arm   = -9;
  short ptbin = -9;
  int blueSpinConfig = -9;
  int yeloSpinConfig = -9;


  _nFills = 0;
  int oldFill=-9;

  int nRuns=0;
  int runlist[300];
  for (int i=0; i<300; i++) runlist[i] = 0;
  ifstream fin(IConsts::RunList);

  while (fin >> runlist[nRuns])  nRuns++;

  fin.close();

  cout << "nRuns = " << nRuns << endl;

  bool runExists = false;

  for (int i=0; i<jpsi->GetEntries(); i++)
    {
      jpsi->GetEntry(i);

      if (i==0) oldFill = FillNumber;

      runExists = false;
      for (int l=0; l<nRuns; l++)
	{
	  if (runlist[l] == RunNumber) runExists = true;
	}
      if (!runExists) continue;


      Arm = diMuPz > 0 ? North : South;

      if (diMuMass < minMs[Arm] || diMuMass > maxMs[Arm]) continue;

      ptbin = -9;
      for (int p=0; p<_NPtBins; p++)
	if (diMuPt > IConsts::PtT[p] && diMuPt < IConsts::PtT[p+1]) ptbin = p;
      if (ptbin < 0) continue;


      blueSpinConfig = BlueSpinConfig(diMuPhi, Arm, SpinB);
      yeloSpinConfig = YeloSpinConfig(diMuPhi, Arm, SpinB);

      _SumPhi   [ptbin][int(blueSpinConfig/2.)][_nFills] += 1.;
      _SumCosPhi[ptbin][int(blueSpinConfig/2.)][_nFills] +=
	fabs(cos(diMuPhi - BlueSpinShift));

      _SumPhi   [ptbin][int(yeloSpinConfig/2.)][_nFills] += 1.;
      _SumCosPhi[ptbin][int(yeloSpinConfig/2.)][_nFills] +=
	fabs(cos(diMuPhi - YeloSpinShift));


      if (diMuCharge == 0)
      	{
      	  _Njpsi[ptbin][blueSpinConfig][_nFills] += 1.;
      	  _Njpsi[ptbin][yeloSpinConfig][_nFills] += 1.;
      	}

      if (diMuCharge != 0)
      	{
      	  _Nbkgr[ptbin][blueSpinConfig][_nFills] += 1.;
      	  _Nbkgr[ptbin][yeloSpinConfig][_nFills] += 1.;
      	}


      if (FillNumber != oldFill)
      	{
      	  _Fill[_nFills++] = oldFill;
      	  oldFill = FillNumber;
      	}
    }

  cout << _nFills << endl;


  _gNjpsi = new TGraph*[_NPtBins][nSPINCONF];

  for (int i=0; i<_NPtBins; i++)
    for (int j=0; j<nSPINCONF; j++)
      _gNjpsi[i][j] = new TGraph(_nFills, _Fill, _Njpsi[i][j]);


  TCanvas **cNjpsi = new TCanvas*[_NPtBins];

  for (int i=0; i<_NPtBins; i++)
    {
      TString title = "Njpsi_";
      title += i;
      cNjpsi[i] = new TCanvas(title.Data(), title.Data(), 1300,700);
      cNjpsi[i]->Divide(4,4);

      for (int j=0; j<nSPINCONF; j++)
      	{
      	  cNjpsi[i]->cd(j+1);
      	  _gNjpsi[i][j]->Draw("AP");
      	}
    }

}
//-----------------------------------------------------------------------------


int IJpsiAN::BlueSpinConfig(double Phi, short Arm, short SpinB)
{
  double aPhi = fabs(Phi - BlueSpinShift);
  short Side  = aPhi < PIby2 ? West : East;
 
  // Here left means y x p > 0, where y is the unit vector along
  // y-axis in PHENIX lab frame, and p is the momentum direction of
  // the beam. Similarly right is y x p < 0. 

  if (Arm == North)  // forward for blue beam
    {
      if (Side == West) // left for blue beam
	{
	  if (SpinB == Up) return Blue_Up_Left_Forward;
	  if (SpinB == Dn) return Blue_Dn_Left_Forward;
	}

      if (Side == East) // right for blue beam
	{
	  if (SpinB == Up) return Blue_Up_Right_Forward;
	  if (SpinB == Dn) return Blue_Dn_Right_Forward;
	}
    }


  if (Arm == South)  // backward for blue beam
    {
      if (Side == West)  // left for blue beam
	{
	  if (SpinB == Up) return Blue_Up_Left_Backward;
	  if (SpinB == Dn) return Blue_Dn_Left_Backward;
	}

      if (Side == East) // right for blue beam
	{
	  if (SpinB == Up) return Blue_Up_Right_Backward;
	  if (SpinB == Dn) return Blue_Dn_Right_Backward;
	}
    }

  cout << "IJpsiAN::CountBlue" << endl;
  cout << "Cannot determin the spin configuration. Crashing ...... " << endl;
  return -1;
}
//-----------------------------------------------------------------------------


int IJpsiAN::YeloSpinConfig(double Phi, short Arm, short SpinY)
{
  double aPhi = fabs(Phi - YeloSpinShift);
  short  Side = aPhi < PIby2 ? West : East;

  // Here left means y x p > 0, where y is the unit vector along
  // y-axis in PHENIX lab frame, and p is the momentum direction of
  // the beam. Similarly right is y x p < 0. 


  if (Arm == South)  // forward for yellow beam
    {
      if (Side == East) // left for yellow beam
	{
	  if (SpinY == Up) return Yelo_Up_Left_Forward;
	  if (SpinY == Dn) return Yelo_Dn_Left_Forward;
	}

      if (Side == West) // right for yellow beam
	{
	  if (SpinY == Up) return Yelo_Up_Right_Forward;
	  if (SpinY == Dn) return Yelo_Dn_Right_Forward;
	}
    }


  if (Arm == North)  // backward for yellow beam
    {
      if (Side == East)  // left for yellow beam
	{
	  if (SpinY == Up) return Yelo_Up_Left_Backward;
	  if (SpinY == Dn) return Yelo_Dn_Left_Backward;
	}

      if (Side == West) // right for yellow beam
	{
	  if (SpinY == Up) return Yelo_Up_Right_Backward;
	  if (SpinY == Dn) return Yelo_Dn_Right_Backward;
	}
    }

  cout << "IJpsiAN::CountYelo" << endl;
  cout << "Cannot determin the spin configuration. Crashing ...... " << endl;
  return -1;
}
//-----------------------------------------------------------------------------


void IJpsiAN::CalculateAccFunc()
{
  for (int i=0; i<_NPtBins; i++)
    for (int j=0; j<nASYMCONF1; j++)
      for (int k=0; k<_nFills; k++)
	{
	  _AccFunc[i][j][k] = _SumPhi[i][j][k]/_SumCosPhi[i][j][k];
	}


  TCanvas *cAcc[2];
  cAcc[0] = new TCanvas("cAcc0","Acceptanc function", 1200, 900);
  cAcc[1] = new TCanvas("cAcc1","Acceptanc function", 1200, 900);

  cAcc[0]->Divide(4,2);
  cAcc[1]->Divide(4,2);

  const char* title1[] = {"Blue_Left_Forward",  "Blue_Right_Forward",
			  "Blue_Left_Backward", "Blue_Right_Backward",
			  "Yelo_Left_Forward",  "Yelo_Right_Forward",
			  "Yelo_Left_Backward", "Yelo_Right_Backward"};

  const char* title2[] = {", 0 < pt < 1.4", "1.4 < pt < 6.0"};

  TGraph *gAcc[_NPtBins][nASYMCONF1];
  for (int i=0; i<_NPtBins; i++)
    for (int j=0; j<nASYMCONF1; j++)
      {
	TString title = title1[j];
	title += title2[i];
	title += ";fill number";

	cAcc[i]->cd(1 + j);
	gAcc[i][j] = new TGraph(_nFills, _Fill, _AccFunc[i][j]);
	gAcc[i][j]->SetMinimum(1.);
	gAcc[i][j]->SetMaximum(2.);
	gAcc[i][j]->SetTitle(title.Data());
	gAcc[i][j]->Draw("AP");
	gAcc[i][j]->Fit("pol0");
      }

}
//-----------------------------------------------------------------------------


void IJpsiAN::CalculateAsym()
{
  cout << "IJpsiAN::CalculateAsym" << endl;

  double Nu=0., Nd=0.;
  double R=0., dR=0., R2=0., dR2=0.;
  double P=0., dP=0., P2=0., dP2=0.;
  double f=0., df=0., f2=0., df2=0.;

  double A1=0., A2=0., wA1=0., wA2=0.;

  int Beam = -9;

  for (int i=0; i<_NPtBins; i++)
    for (int k=0; k<_nFills; k++)
      {
	for (int j=0; j<nASYMCONF1; j++)
	  {
	    Nu = _Njpsi[i][0 + j*2][k];
	    Nd = _Njpsi[i][1 + j*2][k];

	    Beam = j < 4 ? Blue : Yelo;

	    R   = _R[Beam][k];
	    dR  = _dR[Beam][k];
	    R2  = R*R;
	    dR2 = dR*dR;

	    P   = _Pol[Beam][k];
	    dP  = _ePolStat[Beam][k];
	    P2  = P*P;
	    dP2 = dP*dP;

	    f   = _AccFunc[i][j][k];
	    df  = 0.;  // not yet calculated
	    f2  = f*f;
	    df2 = df*df;

	    _AN[i][j][k] = f*(Nu - R*Nd)/(Nu + R*Nd)/P;

	    _dAN[i][j][k] = (1/P2/pow(Nu + R*Nd, 2))*
	      sqrt( (dP2*f2 + df2*P2)*pow(Nu*Nu - Nd*Nd*R2, 2) +
		    4.*P2*f2*Nu*Nd* (Nu*Nd*dR2 + (Nu + Nd)*R2) );

	    cout << i << "  " << j << "  " << k << "  " << _Fill[k] << "  "
		 << _AN[i][j][k] << "   " << _dAN[i][j][k] << endl;
	    cout << "R, P  " << R << "  " << dR << "  " 
		 << P << "  " << dP << endl;
	  }


	// weighted averages.
	for (int j=0; j<4; j++)
	  {
	    A1  = _AN[i][0 + j*2][k];
	    wA1 = 1./
	      (_dAN[i][0 + j*2][k]*_dAN[i][0 + j*2][k]);

	    A2  = _AN[i][1 + j*2][k];
	    wA2 = 1./
	      (_dAN[i][1 + j*2][k]*_dAN[i][1 + j*2][k]);

	    _AN[i][8+j][k]  = (A1 * wA1  - A2 * wA2)/(wA1 + wA2);
	    _dAN[i][8+j][k] = sqrt(1./(wA1 + wA2));
	  }


	A1  = _AN[i][Blue_Forward][k];
	wA1 = 1./
	  (_dAN[i][Blue_Forward][k]*_dAN[i][Blue_Forward][k]);

	A2  = _AN[i][Yelo_Forward][k];
	wA2 = 1./
	  (_dAN[i][Yelo_Forward][k]*_dAN[i][Yelo_Forward][k]);

	_AN[i][Forward][k]  = (A1 * wA1  + A2 * wA2)/(wA1 + wA2);
	_dAN[i][Forward][k] = sqrt(1./(wA1 + wA2));



	A1  = _AN[i][Blue_Backward][k];
	wA1 = 1./
	  (_dAN[i][Blue_Backward][k]*_dAN[i][Blue_Backward][k]);

	A2  = _AN[i][Yelo_Backward][k];
	wA2 = 1./
	  (_dAN[i][Yelo_Backward][k]*_dAN[i][Yelo_Backward][k]);

	_AN[i][Backward][k]  = (A1 * wA1  + A2 * wA2)/(wA1 + wA2);
	_dAN[i][Backward][k] = sqrt(1./(wA1 + wA2));
      }


  TCanvas *cAN[4];
  cAN[0] = new TCanvas("cAN0","A_N0", 1200, 900);
  cAN[1] = new TCanvas("cAN1","A_N1", 1200, 900);
  cAN[2] = new TCanvas("cAN2","A_N2", 1200, 900);
  cAN[3] = new TCanvas("cAN3","A_N3", 1200, 900);

  cAN[0]->Divide(4,2);
  cAN[1]->Divide(4,2);
  cAN[2]->Divide(3,2);
  cAN[3]->Divide(3,2);


  const char* title1[] = {"Blue_Left_Forward",  "Blue_Right_Forward",
  			  "Blue_Left_Backward", "Blue_Right_Backward",
  			  "Yelo_Left_Forward",  "Yelo_Right_Forward",
  			  "Yelo_Left_Backward", "Yelo_Right_Backward",
			  "Blue_Forward", "Blue_Backward",
			  "Yelo_Forward", "Yelo_Backward",
			  "Forward", "Backward"};

  const char* title2[] = {", 0 < pt < 1.4", ", 1.4 < pt < 6.0"};

  TGraphErrors *gAN[_NPtBins][nASYMCONF2];
  for (int i=0; i<_NPtBins; i++)
    for (int j=0; j<nASYMCONF2; j++)
      {
  	TString title = title1[j];
  	title += title2[i];
  	title += ";fill number";

  	cAN[i]->cd(1 + j);
  	gAN[i][j] = new TGraphErrors(_nFills, _Fill, _AN[i][j], 0, _dAN[i][j]);
  	gAN[i][j]->SetMinimum(-4.);
  	gAN[i][j]->SetMaximum( 4.);
  	gAN[i][j]->SetTitle(title.Data());

	if (j<8) cAN[i]->cd(1+j);
	else cAN[i+2]->cd(j-7);

	gAN[i][j]->Draw("AP");
	gAN[i][j]->Fit("pol0");
      }

}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
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
    N1*(5./100)*TMath::Gaus(x, M1, S1, true) +
    N2*(5./100)*TMath::Gaus(x, M2, S2, true) ;

  return f;
}
//-----------------------------------------------------------------------------
