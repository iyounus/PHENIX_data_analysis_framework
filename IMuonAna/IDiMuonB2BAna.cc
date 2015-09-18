#include "IDiMuonB2BAna.hh"
#include <iostream>
#include <fstream>
#include <assert.h>

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "IPdst.hh"
#include "IEvent.hh"
#include "IConsts.hh"

#include "IHeader.hh"
#include "IHeaderList.hh"

#include "IMuon.hh"
#include "IMuonList.hh"

#include "IFillLookup.hh"
#include "ISpinPattern.hh"

using namespace std;

#define MINMS  1.0
#define MAXMS  16.0

IDiMuonB2BAna::IDiMuonB2BAna(const char* file)
{
  assert(IConsts::Defined); // IConsts must be defined before calling this ctor
  IConsts::Print();

  _vEvent = new vector<IEvent *>;
  _vHist = new vector<TH1D *>;

  _outFile = new TFile(file, "recreate");

  _dimu = new TTree("dimuons", "di muon ana");

  _dimu->Branch("RunNumber",      &_RunNumber,      "RunNumber/I");
  _dimu->Branch("FillNumber",     &_FillNumber,     "FillNumber/I");
  _dimu->Branch("EvtNumber",      &_EvtNumber,      "EvtNumber/I");
  _dimu->Branch("NMuons",         &_NMuons,         "NMuons/S");
  _dimu->Branch("SpinXingID",     &_SpinXingID,     "SpinXingID/S");
  _dimu->Branch("SpinB",          &_SpinB,          "SpinB/S");
  _dimu->Branch("SpinY",          &_SpinY,          "SpinY/S");
  _dimu->Branch("Zvtx",           &_Zvtx,           "Zvtx/D");

  _dimu->Branch("Mu1_DDG0",       &_Mu1_DDG0,       "Mu1_DDG0/D");
  _dimu->Branch("Mu1_DG0",        &_Mu1_DG0,        "Mu1_DG0/D");
  _dimu->Branch("Mu1_DS3",        &_Mu1_DS3,        "Mu1_DS3/D");
  _dimu->Branch("Mu1_DS3ctp",     &_Mu1_DS3ctp,     "Mu1_DS3ctp/D");
  _dimu->Branch("Mu1_MuTrChi2",   &_Mu1_MuTrChi2,   "Mu1_MuTrChi2/D");
  _dimu->Branch("Mu1_MuIdChi2",   &_Mu1_MuIdChi2,   "Mu1_MuIdChi2/D");
  _dimu->Branch("Mu1_Px",         &_Mu1_Px,         "Mu1_Px/D");
  _dimu->Branch("Mu1_Py",         &_Mu1_Py,         "Mu1_Py/D");
  _dimu->Branch("Mu1_Pz",         &_Mu1_Pz,         "Mu1_Pz/D");
  _dimu->Branch("Mu1_Pt",         &_Mu1_Pt,         "Mu1_Pt/D");
  _dimu->Branch("Mu1_P",          &_Mu1_P,          "Mu1_P/D");
  _dimu->Branch("Mu1_Charge",     &_Mu1_Charge,     "Mu1_Charge/S");

  _dimu->Branch("Mu1_nMuTrHits",  &_Mu1_nMuTrHits,  "Mu1_nMuTrHits/S");
  _dimu->Branch("Mu1_nMuIdHits",  &_Mu1_nMuIdHits,  "Mu1_nMuIdHits/S");
  _dimu->Branch("Mu1_MuTrHitPat", &_Mu1_MuTrHitPat, "Mu1_MuTrHitPat/i");
  _dimu->Branch("Mu1_MuIdHitPat", &_Mu1_MuIdHitPat, "Mu1_MuIdHitPat/i");

  _dimu->Branch("Mu2_DDG0",       &_Mu2_DDG0,       "Mu2_DDG0/D");
  _dimu->Branch("Mu2_DG0",        &_Mu2_DG0,        "Mu2_DG0/D");
  _dimu->Branch("Mu2_DS3",        &_Mu2_DS3,        "Mu2_DS3/D");
  _dimu->Branch("Mu2_DS3ctp",     &_Mu2_DS3ctp,     "Mu2_DS3ctp/D");
  _dimu->Branch("Mu2_MuTrChi2",   &_Mu2_MuTrChi2,   "Mu2_MuTrChi2/D");
  _dimu->Branch("Mu2_MuIdChi2",   &_Mu2_MuIdChi2,   "Mu2_MuIdChi2/D");
  _dimu->Branch("Mu2_Px",         &_Mu2_Px,         "Mu2_Px/D");
  _dimu->Branch("Mu2_Py",         &_Mu2_Py,         "Mu2_Py/D");
  _dimu->Branch("Mu2_Pz",         &_Mu2_Pz,         "Mu2_Pz/D");
  _dimu->Branch("Mu2_Pt",         &_Mu2_Pt,         "Mu2_Pt/D");
  _dimu->Branch("Mu2_P",          &_Mu2_P,          "Mu2_P/D");
  _dimu->Branch("Mu2_Charge",     &_Mu2_Charge,     "Mu2_Charge/S");

  _dimu->Branch("Mu2_nMuTrHits",  &_Mu2_nMuTrHits,  "Mu2_nMuTrHits/S");
  _dimu->Branch("Mu2_nMuIdHits",  &_Mu2_nMuIdHits,  "Mu2_nMuIdHits/S");
  _dimu->Branch("Mu2_MuTrHitPat", &_Mu2_MuTrHitPat, "Mu2_MuTrHitPat/i");
  _dimu->Branch("Mu2_MuIdHitPat", &_Mu2_MuIdHitPat, "Mu2_MuIdHitPat/i");

  _dimu->Branch("diMuCharge",     &_diMuCharge,     "diMuCharge/S");
  _dimu->Branch("diMuMass",       &_diMuMass,       "diMuMass/D");
  _dimu->Branch("diMuP",          &_diMuP,          "diMuP/D");
  _dimu->Branch("diMuPt",         &_diMuPt,         "diMuPt/D");
  _dimu->Branch("diMuPz",         &_diMuPz,         "diMuPz/D");
  _dimu->Branch("diMuPhi",        &_diMuPhi,        "diMuPhi/D");
  _dimu->Branch("diMuEta",        &_diMuEta,        "diMuEta/D");


  CreateHistos();

  _fillLookup = new IFillLookup();
  _spin = new ISpinPattern(IConsts::SpinFile, IConsts::RunList, "Run");

  ifstream fin(IConsts::RunList);
  int run;

  int nrun=0;
  while (fin >> run)
    {
      if (!_spin->Exists(run))
	{
	  cout << " ****** No spin information exists for run "
	       << run << endl;
	  continue;
	}

      _spin->GetSpinB(run, _spinB);
      _spin->GetSpinY(run, _spinY);

      cout << ++nrun << "\t";
      ReadPdst(run);

      FillHistos();

      ResetEventVector();
    }


  _outFile->Write();
  _outFile->Close();
}
//-----------------------------------------------------------------------------


IDiMuonB2BAna::~IDiMuonB2BAna()
{
  cout << "IDiMuonB2BAna::~IDiMuonB2BAna" << endl;

  for (unsigned int i=0; i<_vEvent->size(); i++)
    delete _vEvent->at(i);

  _vEvent->clear();
  delete _vEvent;

  for (unsigned int i=0; i<_vHist->size(); i++)
    delete _vHist->at(i);

  _vHist->clear();
  delete _vHist;
}
//-----------------------------------------------------------------------------


void IDiMuonB2BAna::CreateHistos()
{
  TString name = "";
  // same order as enum SIGN
  char *sign[] = {(char*)"MM", (char*)"PM", (char*)"PP"};

  for (int i=0; i<3; i++)
    {
      name = "hDiMuMass_";
      name += sign[i];

      _hDiMuMass[i] = new TH1D(name.Data(), name.Data(), 150, MINMS, MAXMS);
      _vHist->push_back(_hDiMuMass[i]);

      name = "hDiMuPt_";
      name += sign[i];

      _hDiMuPt[i] = new TH1D(name.Data(), name.Data(), 100, 0, 20);
      _vHist->push_back(_hDiMuPt[i]);
    }
}
//-----------------------------------------------------------------------------


void IDiMuonB2BAna::ReadPdst(int run)
{
  _pdst = new IPdst(run);
  cout << "Number of entries  " << _pdst->GetEntries() << endl;

  double pz=0., dg0=0.;

  for (int e=0; e<_pdst->GetEntries(); e++)
    {
      _pdst->GetEntry(e);

      if (_pdst->GetNMuons() < 2) continue; // at least 2 muons should exist

      _Zvtx = _pdst->GetHeader()->GetZVertex();

      //if (fabs(_Zvtx) > 35.) continue;

      IEvent *evt = new IEvent();

      for (short i=0; i<_pdst->GetNMuons(); i++)
       	{
       	  // single muon cuts go here
       	  IMuon *mu = _pdst->GetMuon(i);

	  pz  = mu->GetPz(IMuon::VTX);
	  dg0 = mu->GetDG0();

	  // all single muon cuts come here
	  if (fabs(pz) < 1.4 || fabs(pz) > 40.) continue;
	  if (pz > 0. && dg0 > 25.) continue;
	  if (pz < 0. && dg0 > 30.) continue;
	  if (mu->GetP(IMuon::VTX) > 80.) continue;
	  if (mu->GetPt(IMuon::VTX) > 20.) continue;
	  if (mu->GetMuTrChi2() > 30.) continue;
	  if (mu->GetDDG0() > 10.) continue;
	  //if (mu->GetMuIdHitPat() < 64) continue;

       	  // if the muon passes all the cuts, add this to event
       	  evt->AddMuon(mu);
       	}

      if (evt->GetNMuons() > 1) // if there are at least 2 muon, keep the event
      	{
      	  evt->AddHeader(_pdst->GetHeader());
      	  _vEvent->push_back(evt);
      	}
      else
	delete evt;

    }// loop over events

  delete _pdst;
  _pdst = NULL;
}
//-----------------------------------------------------------------------------


void IDiMuonB2BAna::FillHistos()
{
  cout << "IDiMuonB2BAna::FillHistos" << endl;

  TLorentzVector *diMu = new TLorentzVector(0.,0.,0.,0.);

  _Side   = -9;
  _PtBin  = -9;
  _diMuCharge = 0;

  short sign;
  short char1, char2;

  for (unsigned int e=0; e<_vEvent->size(); e++)
    {
      IEvent *evt = _vEvent->at(e);
      _NMuons = evt->GetNMuons();

      IHeader *hdr = evt->GetHeader();

      _EvtNumber = hdr->GetEventID();
      _Zvtx      = hdr->GetZVertex();

      _SpinXingID = hdr->GetSpinGL1CrossingID();
      _SpinB = _spinB[_SpinXingID];
      _SpinY = _spinY[_SpinXingID];

      if (_SpinB == 0 || _SpinY == 0) continue; // remove empty bunches

      _RunNumber  = hdr->GetRunID();
      _FillNumber = _fillLookup->GetFill(_RunNumber);

      for (int i=0; i<_NMuons-1; i++)
	{
	  IMuon *Mu1 = evt->GetMuon(i);

	  char1 = Mu1->GetCharge();

	  for(int j=i+1; j<_NMuons; j++)
	    {
	      IMuon *Mu2 = evt->GetMuon(j);

	      char2 = Mu2->GetCharge();

	      // both muons in opposite arms
	      if (Mu1->GetPz(IMuon::VTX)*Mu2->GetPz(IMuon::VTX) > 0) continue;


	      // 0->both negative, 1->opposite signs, 2-> both positive
	      sign = short(1+(char1+char2)/2);

	      // all dimuon cuts go here

	      diMu->SetPxPyPzE(Mu1->GetPx(IMuon::VTX) + Mu2->GetPx(IMuon::VTX),
			       Mu1->GetPy(IMuon::VTX) + Mu2->GetPy(IMuon::VTX),
			       Mu1->GetPz(IMuon::VTX) + Mu2->GetPz(IMuon::VTX),
			       Mu1->GetP (IMuon::VTX) + Mu2->GetP (IMuon::VTX));

	      _diMuMass   = diMu->M();
	      _diMuEta    = diMu->Rapidity();
	      _diMuPt     = diMu->Pt();

	      if (_diMuMass < MINMS || _diMuMass > MAXMS) continue;
	      //if (fabs(_diMuEta) < 1.2 || fabs(_diMuEta) > 2.2) continue;
	      //if (_diMuPt == 0. || _diMuPt > 6.) continue;

	      _diMuP      = diMu->P();
	      _diMuPz     = diMu->Pz();
	      _diMuPhi    = diMu->Phi();
	      _diMuPt     = diMu->Perp();
	      _diMuCharge = char1 + char2;


	      _hDiMuMass[sign]->Fill(_diMuMass);
	      _hDiMuPt[sign]->Fill(_diMuPt);


	      _Mu1_DDG0       = Mu1->GetDDG0();
	      _Mu1_DG0        = Mu1->GetDG0();
	      _Mu1_DS3        = Mu1->GetDS3();
	      _Mu1_DS3ctp     = Mu1->GetDS3cpt();
	      _Mu1_MuTrChi2   = Mu1->GetMuTrChi2();
	      _Mu1_MuIdChi2   = Mu1->GetMuIdChi2();
	      _Mu1_Px         = Mu1->GetPx(IMuon::VTX);
	      _Mu1_Py         = Mu1->GetPy(IMuon::VTX);
	      _Mu1_Pz         = Mu1->GetPz(IMuon::VTX);
	      _Mu1_Pt         = Mu1->GetPt(IMuon::VTX);
	      _Mu1_P          = Mu1->GetP(IMuon::VTX);
	      _Mu1_Charge     = Mu1->GetCharge();

	      _Mu1_nMuTrHits  = Mu1->GetNMuTrHits();
	      _Mu1_nMuIdHits  = Mu1->GetNMuIdHits();
	      _Mu1_MuTrHitPat = Mu1->GetMuTrHitPat();
	      _Mu1_MuIdHitPat = Mu1->GetMuIdHitPat();


	      _Mu2_DDG0       = Mu2->GetDDG0();
	      _Mu2_DG0        = Mu2->GetDG0();
	      _Mu2_DS3        = Mu2->GetDS3();
	      _Mu2_DS3ctp     = Mu2->GetDS3cpt();
	      _Mu2_MuTrChi2   = Mu2->GetMuTrChi2();
	      _Mu2_MuIdChi2   = Mu2->GetMuIdChi2();
	      _Mu2_Px         = Mu2->GetPx(IMuon::VTX);
	      _Mu2_Py         = Mu2->GetPy(IMuon::VTX);
	      _Mu2_Pz         = Mu2->GetPz(IMuon::VTX);
	      _Mu2_Pt         = Mu2->GetPt(IMuon::VTX);
	      _Mu2_P          = Mu2->GetP(IMuon::VTX);
	      _Mu2_Charge     = Mu2->GetCharge();

	      _Mu2_nMuTrHits  = Mu2->GetNMuTrHits();
	      _Mu2_nMuIdHits  = Mu2->GetNMuIdHits();
	      _Mu2_MuTrHitPat = Mu2->GetMuTrHitPat();
	      _Mu2_MuIdHitPat = Mu2->GetMuIdHitPat();

	      _dimu->Fill();
	    }
	}
    }
}
//-----------------------------------------------------------------------------


void IDiMuonB2BAna::ResetEventVector()
{
  for (unsigned int i=0; i<_vEvent->size(); i++)
    delete _vEvent->at(i);

  _vEvent->clear();
}
//-----------------------------------------------------------------------------
