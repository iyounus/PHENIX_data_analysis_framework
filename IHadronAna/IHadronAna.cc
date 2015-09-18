#include "IHadronAna.hh"

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

#include "ITrack.hh"
#include "ITrackList.hh"

using namespace std;

IHadronAna::IHadronAna()
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

  while (fin >> runNum)
    {
      ReadPdst(runNum);
      FillHistos();
      ResetEventVector();
    }

  _outFile->Write();
  _outFile->Close();
}
//-----------------------------------------------------------------------------


IHadronAna::~IHadronAna()
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


void IHadronAna::ResetEventVector()
{
  cout << "IHadronAna::ResetEventVector" << endl;

  for (unsigned int i=0; i<_vEvent->size(); i++)
    delete _vEvent->at(i);

  _vEvent->clear();
  // NOTE since _vEvent is a pointer to a vector of pointers, the clear() 
  // function does not call the destructor of IEvent class. So, I have to first
  // manually delete every IEvent object and then call clear to reset vector.
}
//-----------------------------------------------------------------------------


void IHadronAna::CreateHistos()
{
  cout << "IHadronAna::CreateHistos" << endl;


  _hPt = new TH1D*[IConsts::NPtTBins];
  _hEP = new TH1D*[IConsts::NPtTBins];

  TString hname = "";
  for (int i=0; i<IConsts::NPtTBins; i++)
    {
      hname = "hPt";
      hname += i;

      _hPt[i] = new TH1D(hname.Data(), hname.Data(),
			 50, IConsts::PtT[i], IConsts::PtT[i+1]);
      _vHist->push_back(_hPt[i]);


      hname = "hEP";
      hname += i;

      _hEP[i] = new TH1D(hname.Data(), hname.Data(), 50, 0., 1.5);
      _vHist->push_back(_hEP[i]);
    }

}
//-----------------------------------------------------------------------------


void IHadronAna::ResetHistos()
{
  for (unsigned int i=0; i<_vHist->size(); i++)
    _vHist->at(i)->Reset();
}
//-----------------------------------------------------------------------------


void IHadronAna::ReadPdst(int run)
{
  _pdst = new IPdst(run);
  cout << "Number of entries  " << _pdst->GetEntries() << endl;

  for (int e=0; e<_pdst->GetEntries(); e++)
    {
      _pdst->GetEntry(e);

      if (_pdst->GetNTracks() < 1) continue;

      IEvent *evt = new IEvent();

      for (short i=0; i<_pdst->GetNTracks(); i++)
       	{
       	  // track cuts go here
       	  ITrack *track = _pdst->GetTrack(i);

	  if (track->GetPt() < IConsts::PtT[0]) continue;
	  if (track->GetPt() > IConsts::PtT[IConsts::NPtTBins]) continue;

	  double pc3s = track->GetPc3MatchingSig();
	  double emcs = track->GetEmcMatchingSig();
	  if (pc3s > IConsts::MatchSig && pc3s < 14000.) continue;
	  if (pc3s > 14000. && emcs > IConsts::MatchSig) continue;

	  short qt = track->GetQuality();
	  if (!(qt == 31 || qt > 60)) continue;

	  if (fabs(track->GetZedDC()) > 75.) continue;
	  if (fabs(track->GetZedDC()) <  5.) continue;


       	  // if the track passes all the cuts, add this to event
       	  evt->AddTrack(track);
       	}


      if (evt->GetNTracks() > 0)
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


void IHadronAna::FillHistos()
{
  cout << "IHadronAna::FillHistos" << endl;

  double pt=0.;

  short ptbin = -9;

  for (unsigned int e=0; e<_vEvent->size(); e++)
    {
      IEvent *evt = _vEvent->at(e);

      for (int i=0; i<evt->GetNTracks(); i++)
	{
	  ITrack *tr = evt->GetTrack(i);

	  pt = tr->GetPt();

	  for (int j=0; j<IConsts::NPtTBins; j++)
	    if (pt > IConsts::PtT[j] && pt < IConsts::PtT[j+1])
	      {
		ptbin = j;
		break;
	      }

	  _hPt[ptbin]->Fill(pt);
	  _hEP[ptbin]->Fill(tr->GetEPratio());
	}
    }
}
//-----------------------------------------------------------------------------
