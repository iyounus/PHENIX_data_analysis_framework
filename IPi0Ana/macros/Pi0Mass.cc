#include "Pi0Mass.hh"

#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TString.h"

using namespace std;

#define NPTBINS 4

Pi0Mass::Pi0Mass(const char *instance, const char* det,
		 const char *file)
{
  TString name = "pi0mass";
  name += instance;
  _hMass = new TH1D(name.Data(),name.Data(), 125, 0, 0.25);

  _f = new TFile(file);
  if (_f->IsZombie()) exit(0);

  //double pt[] = {2.0, 2.5, 3.0, 3.5, 4.2, 5.2, 6.4, 8.0, 10.0};

  for (int i=0; i<NPTBINS; i++)
    {
      name = "pi0Mass_";
      name += det;
      name += i;

      _hmass[i] = (TH1D *)_f->Get(name.Data());

      name = "pt_";
      name += det;
      name += i;

      _hpt[i] = (TH1D *)_f->Get(name.Data());
      _meanPt[i] = _hpt[i]->GetMean();
    }


  for (int i=1; i<NPTBINS; i++) // ignore first 2 pt bins
    _hMass->Add(_hmass[i]);


  double low=0.06, high=0.25;
  for (int i=0; i<NPTBINS; i++)
    {
      name = "pi0peak";
      name += instance;
      name += i;
      _fn[i] = new TF1(name.Data(),"gaus(0) + pol3(3)", low, high);

      name = "signal";
      name += instance;
      name += i;
      _sign[i] = new TF1(name.Data(),"gaus(0)",low, high);

      _sign[i]->SetParNames("Const","Mean","Sigma");

      name = "bkgr";
      name += instance;
      name += i;
      _bkgr[i] = new TF1(name.Data(),"pol3(0)",low, high);

//       Fit_hmass(i, "QRNO");

//       double mean  = _fn[i]->GetParameter(1);
//       double sigma =  fabs(_fn[i]->GetParameter(2));

//       //cout << "\n" << mean << "\t" << sigma << "\t";

//       double min = mean - 2.* sigma;
//       double max = mean + 2.* sigma;

//       _SS[i] = _sign[i]->Integral(min, max);
//       _BB[i] = _bkgr[i]->Integral(min, max);
//       _SB[i] = _SS[i]/_BB[i];

//       //cout << _SS[i] << "\t" << _BB[i] << "\t" << _SB[i] << "\t";

//       _nPi0s[i] = _hmass[i]->Integral(_hmass[i]->FindBin(min),
// 				      _hmass[i]->FindBin(max));

//       cout << _nPi0s[i] << "\t-------------------\n" << endl;

    }

  _gSB = new TGraph(NPTBINS, _meanPt, _SB);
  _gSB->SetTitle("Signal/Bkgr;mean Pt");

}
//-----------------------------------------------------------------------------


Pi0Mass::~Pi0Mass()
{
  cout << "Pi0Mass::~Pi0Mass" << endl;
  _f->Close();
  delete _f;

  for (int i=0; i<NPTBINS; i++)
    {
      delete _fn[i];
      delete _sign[i];
      delete _bkgr[i];
    }


  delete _hMass;

  //delete _gSB;
}
//-----------------------------------------------------------------------------


void Pi0Mass::Draw_hmass(int i)
{
  _hmass[i]->Draw();
  _fn[i]->Draw("same");
  _bkgr[i]->Draw("same");
}
//-----------------------------------------------------------------------------


void Pi0Mass::Fit_hmass(int i, char *opt)
{
  _fn[i]->
    SetParameter(0, _hmass[i]->Integral(_hmass[i]->FindBin(0.12), 
					_hmass[i]->FindBin(0.158),"width"));
  _fn[i]->SetParameter(1, 0.136);
  _fn[i]->SetParameter(2, 0.085);

  _hmass[i]->Fit(_fn[i], opt);

  _sign[i]->FixParameter(0, _fn[i]->GetParameter(0));
  _sign[i]->FixParameter(1, _fn[i]->GetParameter(1));
  _sign[i]->FixParameter(2, _fn[i]->GetParameter(2));


  _bkgr[i]->FixParameter(0, _fn[i]->GetParameter(3));
  _bkgr[i]->FixParameter(1, _fn[i]->GetParameter(4));
  _bkgr[i]->FixParameter(2, _fn[i]->GetParameter(5));
  _bkgr[i]->FixParameter(3, _fn[i]->GetParameter(6));
  _bkgr[i]->SetLineColor(4);
}
//-----------------------------------------------------------------------------


TF1 *Pi0Mass::Fit_hMass(char *opt, double mean1, double sig1)
{
  //cout << "Pi0Mass::Fit_hMass" << endl;

  int i=9;
  double low=0.05, high=0.25;

  _hMass->SetTitle("3 < p_{T}(#pi^{o}) < 10 GeV;di-photon invariant mass [GeV];");

  TString name = "pi0peak";
  name += "_3_10";
  _fn[i] = new TF1(name.Data(),"gaus(0) + pol3(3)", low, high);

  _fn[i]->SetParNames("Const","Mean","Sigma","p0","p1","p2","p3");

  name = "signal";
  name += i+1;
  _sign[i] = new TF1(name.Data(),"gaus(0)",low, high);

  name = "bkgr";
  name += i+1;
  _bkgr[i] = new TF1(name.Data(),"pol3(0)",low, high);


  _fn[i]->
    SetParameter(0, _hMass->Integral(_hMass->FindBin(0.113), 
				     _hMass->FindBin(0.158), "width"));
  _fn[i]->SetParameter(1, mean1);
  _fn[i]->SetParameter(2, sig1);

  _hMass->Fit(_fn[i], opt);

  _sign[i]->FixParameter(0, _fn[i]->GetParameter(0));
  _sign[i]->FixParameter(1, _fn[i]->GetParameter(1));
  _sign[i]->FixParameter(2, _fn[i]->GetParameter(2));
  //_sign[i]->Draw("same");


  _bkgr[i]->FixParameter(0, _fn[i]->GetParameter(3));
  _bkgr[i]->FixParameter(1, _fn[i]->GetParameter(4));
  _bkgr[i]->FixParameter(2, _fn[i]->GetParameter(5));
  _bkgr[i]->FixParameter(3, _fn[i]->GetParameter(6));
  _bkgr[i]->SetLineColor(4);
  _bkgr[i]->Draw("same");


  double mean  = _fn[i]->GetParameter(1);
  double sigma =  fabs(_fn[i]->GetParameter(2));

  //cout << "\n" << mean << "\t" << sigma << "\t";

  double min = mean - 2.* sigma;
  double max = mean + 2.* sigma;

  _SS[i] = _sign[i]->Integral(min, max);
  _BB[i] = _bkgr[i]->Integral(min, max);
  _SB[i] = _SS[i]/_BB[i];

  cout << _SS[i] << "\t" << _BB[i] << "\t" << _SB[i] << endl;

  _nPi0s[i] = _hMass->Integral(_hMass->FindBin(min), _hMass->FindBin(max));

  return _fn[i];
}
//-----------------------------------------------------------------------------
