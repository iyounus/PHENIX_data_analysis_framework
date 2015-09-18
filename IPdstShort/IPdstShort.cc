#include "IPdstShort.hh"

#include <iostream>
#include <fstream>
#include <assert.h>


#include "IConsts.hh"
#include "IPdst.hh"
#include "IHeader.hh"
#include "IPhoton.hh"
#include "ITrack.hh"
#include "IElectron.hh"
#include "IMuon.hh"
#include "IBbc.hh"
#include "IZdcSmd.hh"
#include "IHeaderList.hh"
#include "IPhotonList.hh"
#include "ITrackList.hh"
#include "IElectronList.hh"
#include "IMuonList.hh"
#include "IBbcList.hh"
#include "IZdcSmdList.hh"


#include "TFile.h"
#include "TTree.h"


using namespace std;

typedef unsigned short ushort;

#define NTOWERS        24768 


IPdstShort::IPdstShort(int run, const char* prefix):
  _cMinNPhotons(1),
  _cMinNTracks(0),
  _cMinNMuons(0),
  _cZVtxCentral(30.),
  _cZVtxMuon(55.),
  _cMaxZedDC(75.),
  _cMinZedDC(5.),
  _cMinTrackPt(0.2),
  _cMaxTrackPt(50.),
  _cMinPhotonEn(0.1),
  _cMaxPhotonEn(50.),
  _cMinPhotonProb(0.02),
  _cMinTrackHighPt(0.2),
  _cMinPhotonHighPt(0.4),
  _cMinEmcTof(-10.),
  _cMaxEmcTof( 15.),
  _bPhoton(true),
  _bTrack(true),
  _bElect(true),
  _bMuon(true),
  _bBbc(true),
  _bZdcSmd(true),
  _ApplyWarnMap(false),
  _EMC_badChan(NULL)
{
  cout << "IPdstShort::IPdstShort" << endl;
  assert(IConsts::Defined);  // IConsts must be defined already.


  _pdst = new IPdst(run);
  long nEntries = _pdst->GetEntries();
  int  runNum   = _pdst->GetRunNumberFromPDST();

  cout << "run number in pdst is " << runNum << endl; // sanity check
  cout << "number of entries: " << nEntries << endl;

  TString fileName = IConsts::OutDir;
  fileName += prefix;
  fileName += runNum;
  fileName += ".root";

  cout << "\nOutput file:\n"
       << fileName.Data() << "\n" << endl;


  _lHeader = new IHeaderList();
  _NewFile = new TFile(fileName.Data(), "RECREATE");

  _NewPdst = new TTree("pDST","pDST for whatever analysis");

  // header branch is always created. but it will not be filled if
  // no other particle branch is created.
  _NewPdst->Branch("header", "IHeaderList", &_lHeader);

  // these are the objects that can be written to the output file if
  // the branch is created.
  _lPhoton = NULL;
  _lTrack  = NULL;
  _lElect  = NULL;
  _lMuon   = NULL;
  _lBbc    = NULL;
  _lZdcSmd = NULL;

  cout << "IPdstShort::IPdstShort  done  -------------------------------------"
       << endl;
}
//-----------------------------------------------------------------------------


IPdstShort::~IPdstShort()
{
  cout << "IPdstShort::~IPdstShort" << endl;
  delete _lHeader;

  if (_bPhoton) delete _lPhoton;
  if (_bTrack)  delete _lTrack;
  if (_bElect)  delete _lElect;
  if (_bMuon)   delete _lMuon;
  if (_bBbc)    delete _lBbc;
  if (_bZdcSmd) delete _lZdcSmd;
}
//-----------------------------------------------------------------------------


void IPdstShort::Init()
{
  cout << "IPdstShort::Init  -------------------------------------------------"
       << endl;
    
  if (_bPhoton)
    {
      _lPhoton = new IPhotonList();
      _NewPdst->Branch("photon",   "IPhotonList",   &_lPhoton);
      cout << "Adding photon branch." << endl;
    }
  else
    _pdst->SetBranchStatus("photonList*", 0);


  if (_bTrack)
    {
      _lTrack  = new ITrackList();
      _NewPdst->Branch("track",    "ITrackList",    &_lTrack);
      cout << "Adding track branch." << endl;
    }
  else
    _pdst->SetBranchStatus("trackList*", 0);

  if (_bElect)
    {
      _lElect  = new IElectronList();
      _NewPdst->Branch("electron", "IElectronList", &_lElect);
      cout << "Adding electron branch." << endl;
    }
  else
    _pdst->SetBranchStatus("electronList*", 0);


  if (_bMuon)
    {
      _lMuon   = new IMuonList();
      _NewPdst->Branch("muon",     "IMuonList",     &_lMuon);
      cout << "Adding muon branch." << endl;
    }
  else
    _pdst->SetBranchStatus("muonList*", 0);


  if (_bBbc)
    {
      _lBbc    = new IBbcList();
      _NewPdst->Branch("bbc",      "IBbcList",      &_lBbc);
      cout << "Adding bbc branch." << endl;
    }
  else
    _pdst->SetBranchStatus("bbcList*", 0);


  if (_bZdcSmd)
    {
      _lZdcSmd = new IZdcSmdList();
      _NewPdst->Branch("zdcsmd",   "IZdcSmdList",   &_lZdcSmd);
      cout << "Adding zdcsmd branch." << endl;
    }
  else
    _pdst->SetBranchStatus("zdcsmdList*", 0);



  cout << "IPdstShort::Init  done  -------------------------------------------"
       << endl;
}
//-----------------------------------------------------------------------------


void IPdstShort::Process()
{
  cout << "IPdstShort::Process  ----------------------------------------------"
       << endl;
  Init();

  ushort nPhotons = 0;
  ushort nTracks = 0;
  ushort nMuons = 0;

  for (long evt=0; evt<_pdst->GetEntries(); evt++)
    {
      if (evt%100000==0) cout << evt << endl;

      _pdst->GetEntry(evt);

      if (_bPhoton)
	{
	  if (_pdst->GetNPhotons() == 0 ||
	      _pdst->GetNPhotons() < _cMinNPhotons) continue;
	  nPhotons = FillPhotons();
	  if (nPhotons < _cMinNPhotons) continue;
	}

      if (_bTrack)
	{
	  if (_pdst->GetNTracks() == 0 ||
	      _pdst->GetNTracks()  < _cMinNTracks)  continue;
	  nTracks  = FillTracks();
	  if (nTracks < _cMinNTracks)  continue;
	}

      if (_bMuon)
	{
	  if (_pdst->GetNMuons()   < _cMinNMuons)   continue;
	  nMuons = FillMuons();
	  if (nMuons < _cMinNMuons)   continue;
	}

      if (_bBbc && _pdst->BbcExists())
	{
	  _lBbc->Reset();
	  IBbc *bbc = _lBbc->AddBbc();
	  bbc->Copy(_pdst->GetBbc());
	}

      if (_bZdcSmd && _pdst->ZdcSmdExists())
	{
	  _lZdcSmd->Reset();
	  IZdcSmd *zdc = _lZdcSmd->AddZdcSmd();
	  zdc->Copy(_pdst->GetZdcSmd());
	}

      _lHeader->Reset();
      IHeader *hdr = _lHeader->AddHeader();
      hdr->Copy(_pdst->GetHeader());

      _NewPdst->Fill();
    }

  _NewFile->Write();
  _NewFile->Close();
  cout << "IPdstShort::Process done  -----------------------------------------"
       << endl;
}
//-----------------------------------------------------------------------------


ushort IPdstShort::FillPhotons()
{
  short twID=-9;
  ushort nHiPtPhotons = 0;
  _lPhoton->Reset();

  float t0 = _pdst->GetBbc()->GetBbcT0();

  float en=0., pT=0., tof=0.;
  bool ert = false;

  for (ushort i=0; i<_pdst->GetNPhotons(); i++)
    {
      IPhoton *photon = _pdst->GetPhoton(i);

      if (_ApplyWarnMap)
	{
	  twID = photon->GetTowerID();
	  if (_EMC_badChan[twID]) continue;
	}

      en = photon->GetEnergy();
      if (en < _cMinPhotonEn || en > _cMaxPhotonEn) continue;

      tof = photon->GetTof() - t0;
      if (tof < _cMinEmcTof || tof > _cMaxEmcTof) continue;


      ert = (photon->Fired4x4a() || photon->Fired4x4c());
      pT = photon->GetPt();
      if (ert && pT > _cMinPhotonHighPt) nHiPtPhotons++;

      IPhoton *ph = _lPhoton->AddPhoton();
      ph->Copy(photon); // every object has copy function!
    }
  if (nHiPtPhotons < 1) _lPhoton->Reset();

  return _lPhoton->GetNPhotons();
}
//-----------------------------------------------------------------------------


ushort IPdstShort::FillTracks()
{
  ushort nTracks = 0;
  ushort nHiPtTracks = 0;
  _lTrack->Reset();
  if (_bElect) _lElect->Reset();

  for (ushort i=0; i<_pdst->GetNTracks(); i++)
    {
      ITrack *track = _pdst->GetTrack(i);

      if (fabs(track->GetZedDC()) < _cMinZedDC) continue;
      if (fabs(track->GetZedDC()) > _cMaxZedDC) continue;

      double pc3s = track->GetPc3MatchingSig();
      double emcs = track->GetEmcMatchingSig();
      if (pc3s > IConsts::MatchSig && pc3s < 14000.) continue;
      if (pc3s > 14000. && emcs > IConsts::MatchSig) continue;

      //short qt = track->GetQuality();
      //if (!(qt == 31 || qt > 60)) continue;

      double pT = double(track->GetPt());

      if (pT > _cMaxTrackPt || pT < _cMinTrackPt) continue;
      if (pT > _cMinTrackHighPt) nHiPtTracks++;

      ITrack *tr = _lTrack->AddTrack();
      tr->Copy(track);

      if (_bElect)
	for (ushort j=0; j<_pdst->GetNElectrons(); j++)
	  {
	    IElectron *el = _pdst->GetElectron(j);

	    if (el->GetTrack() == i)
	      {
		IElectron *el2 = _lElect->AddElectron();
		el2->Copy(el);
		el2->SetTrack(nTracks);
	      }
	  }

      nTracks++;
    }
  if (nHiPtTracks == 0)
    {
      _lTrack->Reset();
      if (_bElect) _lElect->Reset();
    }

  return _lTrack->GetNTracks();
}
//-----------------------------------------------------------------------------


ushort IPdstShort::FillMuons()
{
  _lMuon->Reset();

  for (ushort i=0; i<_pdst->GetNMuons(); i++)
    {
      //muon = pdst->GetMuon(i); // not need unless we need to apply some cuts

      IMuon *mu = _lMuon->AddMuon();
      mu->Copy(_pdst->GetMuon(i));
    }
  return _lMuon->GetNMuons();
}
//-----------------------------------------------------------------------------


void IPdstShort::SetEMCalWarnMap(const char* warnmap)
{
  _ApplyWarnMap = true;
  // read bad towers from the fill
  cout << "IPdstReco::ReadEMCalWarnMap" << endl;
  cout << "Reading EMCal warnmap from " << endl;
  cout << warnmap << endl;

  _EMC_badChan = new bool[NTOWERS];
  for (int i=0; i<NTOWERS; i++)
    _EMC_badChan[i] = false;

  ifstream fin(warnmap);
  int towerID = 0, badFlag = 0, nChan = 0;

  while (fin >> towerID >> badFlag)
    if (badFlag > 0)
      {
        _EMC_badChan[towerID] = true;
        nChan++;
      }

  fin.close();

  cout << "Number of masked (hot/dead/edge) EMCal towers: "
       << nChan << endl;
}
//-----------------------------------------------------------------------------
