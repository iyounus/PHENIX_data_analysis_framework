#include "IMuonAna.hh"

#include <iostream>
#include <fstream>
#include <assert.h>

#include "TH1D.h"
#include "TFile.h"
#include "TString.h"

#include "IPdst.hh"
#include "IConsts.hh"
#include "IEvent.hh"

#include "IHeader.hh"
#include "IHeaderList.hh"

#include "IMuon.hh"
#include "IMuonList.hh"

using namespace std;

IMuonAna::IMuonAna()
{
  assert(IConsts::Defined); // IConsts must be defined before calling this ctor
  IConsts::Print();

  _vEvent = new vector<IEvent *>;
  _vHist = new vector<TH1D *>;

  TString outfile = IConsts::OutDir;
  outfile += IConsts::OutFile;

  _outFile = new TFile(outfile.Data(), "recreate");
  CreateHistos();

  ifstream fin(IConsts::RunList);
  int runNum;
  int nruns=0;
  while (fin >> runNum)
    {
      cout << ++nruns << "\t";
      ReadPdst(runNum);
      FillHistos();
      ResetEventVector();
    }

  _outFile->Write();
  _outFile->Close();
}
//-----------------------------------------------------------------------------


IMuonAna::~IMuonAna()
{
  delete _pdst;

  // delete all the events in the memory if any
  ResetEventVector();
  delete _vEvent;

  // histograms are deleted by TFile::Close();
  _vHist->clear();
  delete _vHist;
}
//-----------------------------------------------------------------------------


void IMuonAna::ResetEventVector()
{
  cout << "IMuonAna::ResetEventVector" << endl;

  for (unsigned int i=0; i<_vEvent->size(); i++)
    delete _vEvent->at(i);

  _vEvent->clear();
  // NOTE since _vEvent is a pointer to a vector of pointers, the clear() 
  // function does not call the destructor of IEvent class. So, I have to first
  // manually delete every IEvent object and then call clear to reset vector.
}
//-----------------------------------------------------------------------------


void IMuonAna::CreateHistos()
{
  cout << "IMuonAna::CreateHistos" << endl;

  _hMuIdHits = new TH1D("hMuIdHits","hMuIdHits",9,3,12);
  _vHist->push_back(_hMuIdHits);

  _hMuTrHits = new TH1D("hMuTrHits","hMuTrHits",10,8,18);
  _vHist->push_back(_hMuTrHits);

  const char* mut[] = {"vtx", "st1", "st2", "st3"};
  const char* det[] = {"_North","_South"};

  TString hname = "";
  for (int i=0; i<4; i++)
    for (int j=0; j<2; j++)
      {
	hname = "hPt_";
	hname += mut[i];
	hname += det[j];
	_hPt[j][i] = new TH1D(hname.Data(), hname.Data(), 100, 0., 10.);
	_vHist->push_back(_hPt[j][i]);

	hname = "hPhi_";
	hname += mut[i];
	hname += det[j];
	_hPhi[j][i] = new TH1D(hname.Data(), hname.Data(), 120, -PI, PI);
	_vHist->push_back(_hPhi[j][i]);

	hname = "hPz_";
	hname += mut[i];
	hname += det[j];
	_hPz[j][i] = new TH1D(hname.Data(), hname.Data(), 100, 2., 22.);
	_vHist->push_back(_hPz[j][i]);
      }


  _hDG0      = new TH1D("hDG0",  "hDG0",  80, 0., 80.);
  _vHist->push_back(_hDG0);

  _hDDG0     = new TH1D("hDDG0", "hDDG0", 50, 0., 1.);
  _vHist->push_back(_hDDG0);

  _hDS3      = new TH1D("hDS3",  "hDS3",  80, 0., 80.);
  _vHist->push_back(_hDS3);

  _hDS3cpt   = new TH1D("hDS3cpt","hDS3cpt",50, 0., 50.);
  _vHist->push_back(_hDS3cpt);
}
//-----------------------------------------------------------------------------


void IMuonAna::ResetHistos()
{
  for (unsigned int i=0; i<_vHist->size(); i++)
    _vHist->at(i)->Reset();
}
//-----------------------------------------------------------------------------


void IMuonAna::ReadPdst(int run)
{
  _pdst = new IPdst(run);
  cout << "Number of entries  " << _pdst->GetEntries() << endl;

  for (int e=0; e<_pdst->GetEntries(); e++)
    {
      _pdst->GetEntry(e);

      if (_pdst->GetNMuons() < 1) continue;

      IEvent *evt = new IEvent();

      for (short i=0; i<_pdst->GetNMuons(); i++)
       	{
       	  // single muon cuts go here
       	  IMuon *mu = _pdst->GetMuon(i);

	  if (fabs(mu->GetPx(IMuon::VTX)) > 50.) continue;
	  if (fabs(mu->GetPy(IMuon::VTX)) > 50.) continue;
	  if (fabs(mu->GetPz(IMuon::VTX)) > 20.) continue;

	  if (mu->GetNMuTrHits() < 12) continue;
	  if (mu->GetNMuIdHits() < 4 ) continue;

       	  // if the muon passes all the cuts, add this to event
       	  evt->AddMuon(mu);
       	}


      if (evt->GetNMuons() > 0)
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


void IMuonAna::FillHistos()
{
  cout << "IMuonAna::FillHistos" << endl;

  int det = -1;
  float pz = 0;
  for (unsigned int e=0; e<_vEvent->size(); e++)
    {
      IEvent *evt = _vEvent->at(e);

      for (int i=0; i<evt->GetNMuons(); i++)
	{
	  IMuon *muon = evt->GetMuon(i);

	  if (muon->GetPz(IMuon::VTX) > 0) det = N;
	  else det = S;

	  if (muon->GetPt(IMuon::VTX) < 10.)
	    _hPt[det][IMuon::VTX]->Fill(muon->GetPt(IMuon::VTX));

	  if (muon->GetPt(IMuon::ST1) < 10.)
	    _hPt[det][IMuon::ST1]->Fill(muon->GetPt(IMuon::ST1));

	  if (muon->GetPt(IMuon::ST2) < 10.)
	    _hPt[det][IMuon::ST2]->Fill(muon->GetPt(IMuon::ST2));

	  if (muon->GetPt(IMuon::ST3) < 10.)
	    _hPt[det][IMuon::ST3]->Fill(muon->GetPt(IMuon::ST3));

	  _hPhi[det][IMuon::VTX]->Fill(muon->GetPhi(IMuon::VTX));
	  _hPhi[det][IMuon::ST1]->Fill(muon->GetPhi(IMuon::ST1));
	  _hPhi[det][IMuon::ST2]->Fill(muon->GetPhi(IMuon::ST2));
	  _hPhi[det][IMuon::ST3]->Fill(muon->GetPhi(IMuon::ST3));

	  pz = fabs(muon->GetPz(IMuon::VTX));
	  if (pz > 2. && pz < 22.) _hPz[det][IMuon::VTX]->Fill(pz);

	  pz = fabs(muon->GetPz(IMuon::ST1));
	  if (pz > 2. && pz < 22.) _hPz[det][IMuon::ST1]->Fill(pz);

	  pz = fabs(muon->GetPz(IMuon::ST2));
	  if (pz > 2. && pz < 22.) _hPz[det][IMuon::ST2]->Fill(pz);

	  pz = fabs(muon->GetPz(IMuon::ST3));
	  if (pz > 2. && pz < 22.) _hPz[det][IMuon::ST3]->Fill(pz);


	  _hMuIdHits->Fill(muon->GetNMuIdHits());
	  _hMuTrHits->Fill(muon->GetNMuTrHits());
	  _hDG0     ->Fill(muon->GetDG0());
	  _hDDG0    ->Fill(muon->GetDDG0());
	  _hDS3     ->Fill(muon->GetDS3());
	  _hDS3cpt  ->Fill(muon->GetDS3cpt());
	}
    }
}
//-----------------------------------------------------------------------------
