#include "IPi0Reco.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>


#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1D.h"


#include "IPdst.hh"
#include "IHeader.hh"
#include "IPhoton.hh"
#include "ITrack.hh"
#include "IConsts.hh"


using namespace std;

#define MAXMS 0.25
#define MINMS 0.


IPi0Reco::IPi0Reco(int runNum)
{
  const int nPtBins = IConsts::NPtTBins;

  _vHistos = new vector<TH1D*>;
  CreateHistos();

  _pdst = new IPdst(runNum);
  cout << _pdst->GetEntries() << endl;


  _nEntries = 0;
  double En1=0., Px1=0., Py1=0., Pz1=0., Ph1=0.;
  double En2=0., Px2=0., Py2=0., Pz2=0., Ph2=0.;
  bool   Tr1=false, Tr2=false;
  short  Tw1=0, Tw2=0;
  bool   PbSc=false, PbGl=false;

  // get run number
  _runNum = _pdst->GetRunNumberFromPDST();

  TLorentzVector *pi0 = new TLorentzVector();
  double pi0Pt = 0.;
  double pi0Ms = 0.;
  double pi0Ph = 0;
  double pi0Th = 0;

  IPhoton *pho1;
  IPhoton *pho2;

  for (long evt=0; evt<_pdst->GetEntries(); evt++)
    //for (long evt=0; evt<1000; evt++)
    {
      if (evt%100000 == 0) cout << evt << endl;

      _pdst->GetEntry(evt);

     if (_pdst->GetNPhotons() < 2) continue;

      for (unsigned short j=0; j<_pdst->GetNPhotons()-1; j++)
	{
	  pho1 = _pdst->GetPhoton(j);

	  Tr1 = pho1->FiredErt();
	  Tw1 = pho1->GetTowerID();


	  En1 = double(pho1->GetEnergy());
	  Px1 = pho1->GetPx();
	  Py1 = pho1->GetPy();
	  Pz1 = pho1->GetPz();

	  Ph1 = pho1->GetPhi();
	  if (Ph1 < -PIby2) Ph1 += TwoPI;

	  PbSc = false;
	  PbGl = false;
	  for (unsigned short k=j+1; k<_pdst->GetNPhotons(); k++)
	    {
	      pho2 = _pdst->GetPhoton(k);

	      Tr2 = pho2->FiredErt();
	      Tw2 = pho2->GetTowerID();

	      if (!Tr1 && !Tr2) continue;

	      En2 = double(pho2->GetEnergy());
	      double alpha = fabs( (En2-En1)/(En2+En1) );
	      if (alpha > IConsts::EAsym) continue;

	      Ph2 = pho2->GetPhi();
	      if (Ph2 < -PIby2) Ph2 += TwoPI;

	      if (Ph1 < PIby2 && Ph2 > PIby2) continue;
	      if (Ph1 > PIby2 && Ph2 < PIby2) continue;

	      if (Tw1 < 15553 && Tw2 < 15553) PbSc = true;
	      if (Tw1 > 15552 && Tw2 > 15552) PbGl = true;

	      Px2 = pho2->GetPx();
	      Py2 = pho2->GetPy();
	      Pz2 = pho2->GetPz();

	      pi0->SetPxPyPzE(Px1+Px2, Py1+Py2, Pz1+Pz2, En1+En2);

	      pi0Pt = pi0->Pt();
	      pi0Ms = pi0->M();
	      pi0Ph = pi0->Phi();
	      pi0Th = pi0->Theta();

	      if (pi0Ms  > MAXMS) continue;
	      if (pi0Ms  < MINMS) continue;
	      if (pi0Pt < IConsts::PtT[0]) continue;
	      if (pi0Pt > IConsts::PtT[nPtBins]) continue;
	      if (pi0Ph < -PIby2) pi0Ph += TwoPI;


	      // determine pt bin for the pion
	      int ptbin = -99;
	      for (int h=0; h<nPtBins; h++)
		if (pi0Pt > IConsts::PtT[h] && pi0Pt < IConsts::PtT[h+1])
		  ptbin = h;

	      if (ptbin < 0) continue;

	      if (PbSc)
		{
		  _hmass_PbSc [ptbin]->Fill(pi0Ms);
		  _hpt_PbSc   [ptbin]->Fill(pi0Pt);
		  _hphi_PbSc  [ptbin]->Fill(pi0Ph);
		  _htheta_PbSc[ptbin]->Fill(pi0Th);
		}

	      if (PbGl)
		{
		  _hmass_PbGl [ptbin]->Fill(pi0Ms);
		  _hpt_PbGl   [ptbin]->Fill(pi0Pt);
		  _hphi_PbGl  [ptbin]->Fill(pi0Ph);
		  _htheta_PbGl[ptbin]->Fill(pi0Th);
		}

	      if (!PbSc && !PbGl)
		{
		  _hmass_ScGl [ptbin]->Fill(pi0Ms);
		  _hpt_ScGl   [ptbin]->Fill(pi0Pt);
		  _hphi_ScGl  [ptbin]->Fill(pi0Ph);
		  _htheta_ScGl[ptbin]->Fill(pi0Th);
		}

	      _nEntries++;
	    } // photon loop
	} // photon loop
    } // event loop

  cout << "=======" << endl;
  cout << _nEntries << endl;


  delete _pdst;

  Write();
}
//-----------------------------------------------------------------------------


IPi0Reco::~IPi0Reco()
{
  if (_WarnMap) delete _WarnMap;
}
//-----------------------------------------------------------------------------


void IPi0Reco::CreateHistos()
{
  const int nPtBins = IConsts::NPtTBins;

  _hmass_PbSc  = new TH1D*[nPtBins];
  _hpt_PbSc    = new TH1D*[nPtBins];
  _hphi_PbSc   = new TH1D*[nPtBins];
  _htheta_PbSc = new TH1D*[nPtBins];

  _hmass_PbGl  = new TH1D*[nPtBins];
  _hpt_PbGl    = new TH1D*[nPtBins];
  _hphi_PbGl   = new TH1D*[nPtBins];
  _htheta_PbGl = new TH1D*[nPtBins];

  _hmass_ScGl  = new TH1D*[nPtBins];
  _hpt_ScGl    = new TH1D*[nPtBins];
  _hphi_ScGl   = new TH1D*[nPtBins];
  _htheta_ScGl = new TH1D*[nPtBins];


  TString name = "";
  TString title = "";
  for (int i=0; i<nPtBins; i++)
    {
      double minPt = IConsts::PtT[i];
      double maxPt = IConsts::PtT[i+1];

      TString ptRange = " ";
      ptRange += int(minPt);
      ptRange += ".";
      ptRange += int(10*minPt - 10*int(minPt));
      ptRange += " < pT < ";
      ptRange += int(maxPt);
      ptRange += ".";
      ptRange += int(10*maxPt - 10*int(maxPt));


      // PbSc ----------------------------------
      name = "pi0Mass_PbSc";
      name += i;
      title = "pi0Mass PbSc";
      title += ptRange;
      _hmass_PbSc[i] = new TH1D(name.Data(), title.Data(), 125, MINMS, MAXMS);
      _vHistos->push_back(_hmass_PbSc[i]);

      name = "pt_PbSc";
      name += i;
      _hpt_PbSc[i] = new TH1D(name.Data(), ptRange.Data(), 50, minPt, maxPt);
      _vHistos->push_back(_hpt_PbSc[i]);

      name = "phi_PbSc";
      name += i;
      _hphi_PbSc[i] = new
	TH1D(name.Data(), ptRange.Data(), 120, -PIby2, 3*PIby2);
      _vHistos->push_back(_hphi_PbSc[i]);

      name = "theta_PbSc";
      name += i;
      _htheta_PbSc[i] = new
	TH1D(name.Data(), ptRange.Data(), 60, 1.0472, 2.0944);
      _vHistos->push_back(_htheta_PbSc[i]);


      // PbGl ----------------------------------
      name = "pi0Mass_PbGl";
      name += i;
      title = "pi0Mass PbGl";
      title += ptRange;
      _hmass_PbGl[i] = new TH1D(name.Data(), title.Data(), 125, MINMS, MAXMS);
      _vHistos->push_back(_hmass_PbGl[i]);

      name = "pt_PbGl";
      name += i;
      _hpt_PbGl[i] = new TH1D(name.Data(), ptRange.Data(), 50, minPt, maxPt);
      _vHistos->push_back(_hpt_PbGl[i]);

      name = "phi_PbGl";
      name += i;
      _hphi_PbGl[i] = new
	TH1D(name.Data(), ptRange.Data(), 120, -PIby2, 3*PIby2);
      _vHistos->push_back(_hphi_PbGl[i]);

      name = "theta_PbGl";
      name += i;
      _htheta_PbGl[i] = new
	TH1D(name.Data(), ptRange.Data(), 60, 1.0472, 2.0944);
      _vHistos->push_back(_htheta_PbGl[i]);


      // ScGl ----------------------------------
      name = "pi0Mass_ScGl";
      name += i;
      title = "pi0Mass ScGl";
      title += ptRange;
      _hmass_ScGl[i] = new TH1D(name.Data(), title.Data(), 125, MINMS, MAXMS);
      _vHistos->push_back(_hmass_ScGl[i]);

      name = "pt_ScGl";
      name += i;
      _hpt_ScGl[i] = new TH1D(name.Data(), ptRange.Data(), 50, minPt, maxPt);
      _vHistos->push_back(_hpt_ScGl[i]);

      name = "phi_ScGl";
      name += i;
      _hphi_ScGl[i] = new
	TH1D(name.Data(), ptRange.Data(), 120, -PIby2, 3*PIby2);
      _vHistos->push_back(_hphi_ScGl[i]);

      name = "theta_ScGl";
      name += i;
      _htheta_ScGl[i] = new
	TH1D(name.Data(), ptRange.Data(), 60, 1.0472, 2.0944);
      _vHistos->push_back(_htheta_ScGl[i]);
    }
}
//-----------------------------------------------------------------------------


void IPi0Reco::Write()
{
  TString outFile = IConsts::OutDir; 
  outFile += "/IPi0Mass_";
  outFile += _runNum;
  outFile += ".root";

  TFile *ff = new TFile(outFile.Data(),"recreate");

  for (unsigned int i=0; i<_vHistos->size(); i++)
    {
      _vHistos->at(i)->Write();
      delete _vHistos->at(i);
    }

  ff->Close();
  _vHistos->clear();
}
//-----------------------------------------------------------------------------


void IPi0Reco::ReadEMCalWarnMap(const char* file)
{
  _WarnMap = new short[24768];
  for (int i=0; i<24768; i++)
    _WarnMap[i] = 0;

  ifstream fin(file);
  int id;
  short dummy;
  short bad=0;

  while(fin >> id >> dummy)
    {
      _WarnMap[id] = dummy;
      if (_WarnMap[id] == 1)
	bad++;
    }
  cout << "No. of masked EMCal channels:  " << bad << endl;
  fin.close();
}
//-----------------------------------------------------------------------------
//=============================================================================
