#include "IGammaHadron.hh"

#include <iostream>
#include <fstream>
#include <assert.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "IPdst.hh"
#include "IEvent.hh"
#include "IHeader.hh"
#include "IPhoton.hh"
#include "ITrack.hh"
#include "IConsts.hh"

using namespace std;

IGammaHadron::IGammaHadron()
{
  cout << "IGammaHadron::IGammaHadron" << endl;

  assert(IConsts::Defined);
  IConsts::Print();

  _nPtT = IConsts::NPtTBins;
  _nPtA = IConsts::NPtABins;

  _vEvent = new vector<IEvent *>;

  // vHist is used only to write all the histos and delete at the end.
  // It doesn't matter what order the hists are pushed in vHist.
  _vHist = new vector<TH1D *>;

  TString outFile = IConsts::OutDir;
  outFile += IConsts::OutFile;

  _f = new TFile(outFile.Data(), "recreate");

  CreateTree();
  CreateHistos();
  ReadFiles();

  _f->Write();
  _f->Close();
}
//------------------------------------------------------------------------------


IGammaHadron::IGammaHadron(int run)
{
  cout << "IGammaHadron::IGammaHadron" << endl;

  assert(IConsts::Defined);
  IConsts::Print();

  _nPtT = IConsts::NPtTBins;
  _nPtA = IConsts::NPtABins;

  _vEvent = new vector<IEvent *>;

  _vHist = new vector<TH1D *>;

  TString outFile = IConsts::OutDir;
  outFile += IConsts::OutFile;

  _f = new TFile(outFile.Data(), "recreate");

  CreateTree();
  CreateHistos();
  ReadFile(run);

  _f->Write();
  _f->Close();
}
//------------------------------------------------------------------------------


IGammaHadron::~IGammaHadron()
{
  ResetVector();
  delete _vEvent;

  for (unsigned int i = 0; i < _vHist->size(); i++)
    delete _vHist->at(i);

  delete _vHist;

  delete _f;
}
//------------------------------------------------------------------------------


void IGammaHadron::CreateTree()
{
  _gam = new TTree("gam","gamma hadron tree");

  _gam->Branch("PhotonPt",  &_phoPt,   "PhotonPt/D");
  _gam->Branch("HadronPt",  &_hadPt,   "HadronPt/D");
  _gam->Branch("PhotonPhi", &_phoPhi,  "PhotonPhi/D");
  _gam->Branch("HadronPhi", &_hadPhi,  "HadronPhi/D");
  _gam->Branch("PhotonEn",  &_phoEn,   "PhotonEn/D");
  _gam->Branch("EP",        &_EP,      "EP/D");
  _gam->Branch("za",        &_za,      "za/D");
  _gam->Branch("xE",        &_xE,      "xE/D");
  _gam->Branch("delPhi",    &_dPhi,    "dPhi/D");
  _gam->Branch("ECone",     &_eCone,   "ECone/D");
  _gam->Branch("HadronCh",  &_hadCh,   "HadronCh/D");
  _gam->Branch("DiHadMs",   &_diHadMs, "DiHadMs/D");
  _gam->Branch("DiPhoMs",   &_diPhoMs, "DiPhoMs/D");
  _gam->Branch("Prob",      &_prob,    "prob/D");
}
//------------------------------------------------------------------------------


void IGammaHadron::CreateHistos()
{
  cout << "IGammaHadron::CreateHistos" << endl;

  _hZt = new TH1D**[_nPtT];
  _hPtT = new TH1D**[_nPtT];
  _hDphiPosH = new TH1D**[_nPtT];
  _hDphiNegH = new TH1D**[_nPtT];

  for (int i = 0; i < _nPtT; i++)
    {
      _hZt[i] = new TH1D*[_nPtA];
      _hPtT[i] = new TH1D*[_nPtA];
      _hDphiPosH[i] = new TH1D*[_nPtA];
      _hDphiNegH[i] = new TH1D*[_nPtA];
    }

  _hEP = new TH1D*[_nPtA];
  _hPtA = new TH1D*[_nPtA];
  _hProb = new TH1D*[_nPtA];

  _hEcone = new TH1D*[_nPtT];
  _hDiPhoMs = new TH1D*[_nPtT];
  _hDiHadMs = new TH1D*[_nPtT];

  TString hName = "";
  TString hTitle1 = "";
  TString hTitle2 = "";

  for (int i = 0; i < IConsts::NPtTBins; i++)
    {
      hTitle1 = "";
      hTitle1 += int(IConsts::PtT[i]);
      hTitle1 += " < ptT < ";
      hTitle1 += int(IConsts::PtT[i + 1]);
      hTitle1 += ",  ";

      for (int j = 0; j < IConsts::NPtABins; j++)
	{
	  hTitle2 = hTitle1.Data();
	  hTitle2 += int(IConsts::PtA[j]);
	  hTitle2 += " < ptA < ";
	  hTitle2 += int(IConsts::PtA[j + 1]);

	  hName = "hZt_";
	  hName += i;
	  hName += "_";
	  hName += j;

	  _hZt[i][j] = new TH1D(hName.Data(), hTitle2.Data(), 100, 0., 5.);
	  _vHist->push_back(_hZt[i][j]);

	  hName = "hPtT_";
	  hName += i;
	  hName += "_";
	  hName += j;

	  _hPtT[i][j] = new TH1D(hName.Data(), hTitle2.Data(), 30,
				 IConsts::PtT[i], IConsts::PtT[i + 1]);
	  _vHist->push_back(_hPtT[i][j]);

	  hName = "hDphiPosH_";
	  hName += i;
	  hName += "_";
	  hName += j;

	  _hDphiPosH[i][j] = new TH1D(hName.Data(), hTitle2.Data(), 120,
				      -PIby2, 3* PIby2 );
	  _vHist->push_back(_hDphiPosH[i][j]);

	  hName = "hDphiNegH_";
	  hName += i;
	  hName += "_";
	  hName += j;

	  _hDphiNegH[i][j] = new TH1D(hName.Data(), hTitle2.Data(), 120,
				      -PIby2, 3* PIby2 );
	  _vHist->push_back(_hDphiNegH[i][j]);
	}

      hName = "hEcone_";
      hName += i;
      _hEcone[i] = new TH1D(hName.Data(), hTitle1.Data(), 50, 0., 5.);
      _vHist->push_back(_hEcone[i]);

      hName = "hDiPhoMs_";
      hName += i;
      _hDiPhoMs[i] = new TH1D(hName.Data(), hTitle1.Data(), 40, 0., 0.4);
      _vHist->push_back(_hDiPhoMs[i]);

      hName = "hDiHadMs_";
      hName += i;
      _hDiHadMs[i] = new TH1D(hName.Data(), hTitle1.Data(), 40, 0., 4.0);
      _vHist->push_back(_hDiHadMs[i]);
    }

  for (int i = 0; i < IConsts::NPtABins; i++)
    {
      hTitle2 = "";
      hTitle2 += int(IConsts::PtA[i]);
      hTitle2 += " < ptA < ";
      hTitle2 += int(IConsts::PtA[i + 1]);

      hName = "hPtA_";
      hName += i;

      _hPtA[i] = new TH1D(hName.Data(), hTitle2.Data(), 50, IConsts::PtA[i],
			  IConsts::PtA[i + 1]);
      _vHist->push_back(_hPtA[i]);

      hName = "hEP_";
      hName += i;

      _hEP[i] = new TH1D(hName.Data(), hTitle2.Data(), 60, 0., 1.5);
      _vHist->push_back(_hEP[i]);

      hName = "hProb_";
      hName += i;

      _hProb[i] = new TH1D(hName.Data(), hTitle2.Data(), 50, 0., 1.);
      _vHist->push_back(_hProb[i]);
    }

  _hPtT_all = new TH1D("hPtT_all", "hPtT_all", 100, IConsts::PtT[0],
		       IConsts::PtT[IConsts::NPtTBins]);
  _hPtA_all = new TH1D("hPtA_all", "hPtA_all", 100, IConsts::PtA[0],
		       IConsts::PtA[IConsts::NPtABins]);

  _vHist->push_back(_hPtT_all);
  _vHist->push_back(_hPtA_all);
}
//------------------------------------------------------------------------------


void IGammaHadron::FillHistos()
{
  cout << "IGammaHadron::FillHistos" << endl;

  int pttbin = -9;
  int ptabin = -9;

  IEvent  *event;
  IPhoton *photon;
  ITrack  *track;

  TVector3 vPhoton(0., 0., 0.);
  TVector3 vTrack(0., 0., 0.);

  for (unsigned int evt=0; evt<_vEvent->size(); evt++)
    {
      event = _vEvent->at(evt);

      for (int i=0; i<event->GetNPhotons(); i++)
	{
	  photon = event->GetPhoton(i);

	  _phoPt = photon->GetPt();
	  _phoEn = photon->GetEnergy();

	  pttbin = -9;
	  for (int k = 0; k < _nPtT; k++)
	    if (_phoPt > IConsts::PtT[k] && _phoPt < IConsts::PtT[k + 1])
	      {
		pttbin = k;
		break;
	      }
	  if (pttbin < 0) continue;

	  _hPtT_all->Fill(_phoPt);

	  // calculate cone energy
	  _eCone = event->ECone(i);
	  if (_eCone > 0.1*photon->GetEnergy()) continue;

	  _hEcone[pttbin]->Fill(_eCone);

	  // see if the prompt photon candidate makes a pi0
	  if (i<event->GetNPhotons()-1)
	    for (int j=i; j<event->GetNPhotons(); j++)
	      {
		_diPhoMs = event->DiPhotonMass(i, j);
		_hDiPhoMs[pttbin]->Fill(_diPhoMs);
	      }

	  _phoPhi = photon->GetPhi();
	  if (_phoPhi < -PIby2) _phoPhi += TwoPI;

	  vPhoton.SetPtThetaPhi(_phoPt, photon->GetTheta(), _phoPhi);

	  //event->_promptPhoton[i] = true; // tag this photon as prompt

	  // now combine this photon with tracks
	  for (int j=0; j<event->GetNTracks(); j++)
	    {
	      track = event->GetTrack(j);

	      _hadPt = track->GetPt();

	      ptabin = -9;
	      for (int k = 0; k < _nPtA; k++)
		if (_hadPt > IConsts::PtA[k] && _hadPt < IConsts::PtA[k + 1])
		  {
		    ptabin = k;
		    break;
		  }
	      if (ptabin < 0) continue;

	      _hadPhi = track->GetPhi0();

	      vTrack.SetPtThetaPhi(_hadPt, track->GetTheta0(), _hadPhi);

	      _dPhi = vPhoton.DeltaPhi(vTrack);
	      if (_dPhi < -PIby2) _dPhi += TwoPI;

	      _prob = track->GetProb();

	      //if (_dPhi < PIby2) continue;
	      //if (_prob > 0.2) continue;
	      //if (phophi < PIby2 && trkphi < PIby2) continue;
	      //if (phophi > PIby2 && trkphi > PIby2) continue;
	      if (_prob < 0) continue;

	      //event->_associatedTrack[j] = true;

	      _hadCh = track->GetCharge();
	      if (_hadCh > 0)
		_hDphiPosH[pttbin][ptabin]->Fill(_dPhi);
	      if (_hadCh < 0)
		_hDphiNegH[pttbin][ptabin]->Fill(_dPhi);


	      _za = _hadPt/_phoPt;
	      _hZt[pttbin][ptabin]->Fill(_za);
	      _hPtT[pttbin][ptabin]->Fill(_phoPt);

	      _xE =  - _za*cos(_dPhi);
	      _EP = track->GetEmcE()/track->GetP();

	      _gam->Fill();
	    }//loop over tracks
	}//loop over photons


      double pta = 0.;
      double mass = 0.;
      for (int i=0; i<event->GetNTracks(); i++)
	{
	  //if (!event->_associatedTrack[i]) continue;

	  track = event->GetTrack(i);

	  pta = track->GetPt();

	  ptabin = -9;
	  for (int k = 0; k < _nPtA; k++)
	    if (pta > IConsts::PtA[k] && pta < IConsts::PtA[k + 1])
	      {
		ptabin = k;
		break;
	      }
	  if (ptabin < 0) continue;

	  _hPtA_all->Fill(pta);
	  _hPtA[ptabin]->Fill(pta);
	  _hEP[ptabin]->Fill(track->GetEmcE() / track->GetP());
	  _hProb[ptabin]->Fill(track->GetProb());

	  if (i<event->GetNTracks() - 1)
	    for (int j=i; j<event->GetNTracks(); j++)
	      {
		mass = event->DiHadronMass(i, j);
		_hDiHadMs[ptabin]->Fill(mass);
	      }
	}
    }
}
//------------------------------------------------------------------------------


void IGammaHadron::ReadFiles()
{
  cout << "IGammaHadron::ReadFiles" << endl;

  cout << IConsts::RunList << endl;
  ifstream fin(IConsts::RunList);

  int runNum = 0;
  int nFiles = 0;
  while (fin >> runNum)
    {
      cout << nFiles++ << "   ";
      ReadFile(runNum);
    }
}
//------------------------------------------------------------------------------


void IGammaHadron::ReadFile(int run)
{
  cout << "IGammaHadron::ReadFile" << endl;

  _pdst = new IPdst(run);
  long nEntries = _pdst->GetEntries();

  int runNum = -9999;
  for (long evt = 0; evt < nEntries; evt++)
    {
      _pdst->GetEntry(evt);

      if (runNum != _pdst->GetHeader()->GetRunID())
	{
	  runNum = _pdst->GetHeader()->GetRunID();
	  cout << "run number in pDST  =  " << runNum << endl;
	  break;
	}
    }

  cout << "no. of entries      =  " << nEntries << endl;

  short Tw = 0, Tg = 0;
  IPhoton *photon = 0;
  ITrack *track = 0;

  short nPhotons = 0, nHiPtPhotons = 0, nTriggers = 0;
  short nTracks = 0, nHiPtTracks = 0;
  short *iPhoton = new short[200];
  short *iTrack = new short[200];

  for (long evt = 0; evt < nEntries; evt++)
    {
      if (evt % 10000 == 0) cout << evt << endl;

      _pdst->GetEntry(evt);

      nPhotons = 0;
      nTriggers = 0;
      nHiPtPhotons = 0;
      for (unsigned int i = 0; i < _pdst->GetNPhotons(); i++)
	{
	  photon = _pdst->GetPhoton(i);

	  Tg = photon->FiredErt();
	  Tw = photon->GetTowerID();
	  if (_pdst->CompareTowerID(Tw)) continue;

	  nTriggers += Tg;

	  double pT = photon->GetPt();
	  if (pT > IConsts::PtT[0]) nHiPtPhotons++;

	  iPhoton[nPhotons++] = i;
	}
      if (nHiPtPhotons < 1) continue;
      if (nTriggers < 1) continue; // some photons must trigger the event

      nTracks = 0;
      nHiPtTracks = 0;
      for (unsigned int i = 0; i < _pdst->GetNTracks(); i++)
	{
	  track = _pdst->GetTrack(i);

	  if (fabs(track->GetZedDC()) <  5.) continue;
	  if (fabs(track->GetZedDC()) > 75.) continue;

	  double pc3s = track->GetPc3MatchingSig();
	  double emcs = track->GetEmcMatchingSig();
	  if (pc3s > IConsts::MatchSig && pc3s < 14000.) continue;
	  if (pc3s > 14000. && emcs > IConsts::MatchSig) continue;

	  short qt = track->GetQuality();
	  if (!(qt == 31 || qt > 60)) continue;

	  double pT = double(track->GetPt());
	  if (pT > IConsts::PtA[IConsts::NPtABins]) continue;
	  if (pT > IConsts::PtA[0]) nHiPtTracks++;

	  iTrack[nTracks++] = i;
	}
      if (nHiPtTracks < 1) continue;


      IEvent *event = new IEvent();

      event->AddHeader(_pdst->GetHeader());

      for (int i = 0; i < nPhotons; i++)
	event->AddPhoton(_pdst->GetPhoton(iPhoton[i]));

      for (int i = 0; i < nTracks; i++)
	event->AddTrack(_pdst->GetTrack(iTrack[i]));

      _vEvent->push_back(event);
    }

  delete _pdst;

  FillHistos();
  ResetVector();

  delete[] iPhoton;
  delete[] iTrack;

  return;
}
//------------------------------------------------------------------------------


void IGammaHadron::ResetVector()
{
  for (unsigned int i = 0; i < _vEvent->size(); i++)
    delete _vEvent->at(i);

  _vEvent->clear();
}
//------------------------------------------------------------------------------
