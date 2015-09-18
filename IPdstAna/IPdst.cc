#include "IPdst.hh"

#include <iostream>
#include <fstream>
#include <assert.h>

#include "TChain.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"


#include "IConsts.hh"
#include "IHeader.hh"
#include "IHeaderList.hh"
#include "IPhoton.hh"
#include "IPhotonList.hh"
#include "ITrack.hh"
#include "ITrackList.hh"
#include "IMuon.hh"
#include "IMuonList.hh"
#include "IDiMuon.hh"
#include "IDiMuonList.hh"
#include "IElectron.hh"
#include "IElectronList.hh"


#include "IBbc.hh"
#include "IBbcList.hh"
#include "IZdcSmd.hh"
#include "IZdcSmdList.hh"

#include "IFillLookup.hh"


IPdst::IPdst(int runNum):
  _nPhotons(0),
  _nTracks(0),
  _nElectrons(0),
  _nMuons(0),
  _nDiMuons(0),
  _nNeutrons(0)
{
  //std::cout << "IPdst::IPdst(int runNum)" << std::endl;
  assert(IConsts::Defined);

  TString file = IConsts::InDir;
  file += "/";
  file += IConsts::FilePrefix;
  file += runNum;
  file += "*root";

  std::cout << "Reading  " << file.Data() << std::endl;

  _pDST = new TChain("pDST");
  _pDST->Add(file.Data());

  SetBranches();
}
//-----------------------------------------------------------------------------


IPdst::IPdst(short n, int *runNum):
  _nPhotons(0),
  _nTracks(0),
  _nElectrons(0),
  _nMuons(0),
  _nDiMuons(0),
  _nNeutrons(0)
{
  //std::cout << "IPdst::IPdst(short n, int *runNum)" << std::endl;
  assert(IConsts::Defined);

  _pDST = new TChain("pDST");

  TString file = "";
  for (short i=0; i<n; i++)
    { 
      file = IConsts::InDir;
      file += "/";
      file += IConsts::FilePrefix;
      file += runNum[i];
      file += "*root";

      std::cout << "Reading  " << file.Data() << std::endl;

      _pDST->Add(file.Data());
    }

  SetBranches();
}
//-----------------------------------------------------------------------------


IPdst::IPdst(const char* runList):
  _nPhotons(0),
  _nTracks(0),
  _nElectrons(0),
  _nMuons(0),
  _nDiMuons(0),
  _nNeutrons(0)
{
  //std::cout << "IPdst::IPdst(const char* runList)" << std::endl;
  assert(IConsts::Defined);

  _pDST = new TChain("pDST");

  ifstream fin(runList);

  int runNum;
  TString file = "";
  while (fin >> runNum)
    {
      file = IConsts::InDir;
      file += "/";
      file += IConsts::FilePrefix;
      file += runNum;
      file += "*root";

      std::cout << "Reading  " << file.Data() << std::endl;

      _pDST->Add(file.Data());
    }

  SetBranches();
}
//-----------------------------------------------------------------------------


void IPdst::SetBranches()
{
  _header     = new IHeader();
  _headerList = new IHeaderList();
  _pDST->SetBranchAddress("header", &_headerList);


  // if photon object exists then add branch
  if (_pDST->GetMaximum("nPhotons") > 0)
    {
      _photon     = new IPhoton*[int(_pDST->GetMaximum("nPhotons"))];
      _photonList = new IPhotonList();

      _pDST->SetBranchAddress("photon", &_photonList);
    }
  else
    {
      _photon     = NULL;
      _photonList = NULL;
    }


  // if track object exists then add branch
  if (_pDST->GetMaximum("nTracks") > 0)
    {
      _track      = new ITrack*[int(_pDST->GetMaximum("nTracks"))];
      _trackList  = new ITrackList();

      _pDST->SetBranchAddress("track",  &_trackList);
    }
  else
    {
      _track     = NULL;
      _trackList = NULL;
    }


  // if electron object exists then add branch
  if (_pDST->GetMaximum("nElectrons") > 0)
    {
      _elect      = new IElectron*[int(_pDST->GetMaximum("nElectrons"))];
      _electList  = new IElectronList();

      _pDST->SetBranchAddress("electron",  &_electList);
    }
  else
    {
      _elect     = NULL;
      _electList = NULL;
    }


  // if muon object exists then add branch
  if (_pDST->GetMaximum("nMuons") > 0)
    {
      _muon      = new IMuon*[int(_pDST->GetMaximum("nMuons"))];
      _muonList  = new IMuonList();

      _pDST->SetBranchAddress("muon",  &_muonList);
    }
  else
    {
      _muon     = NULL;
      _muonList = NULL;
    }


  // if dimuon object exists then add branch
  if (_pDST->GetMaximum("nDiMuons") > 0)
    {
      _dimuon      = new IDiMuon*[int(_pDST->GetMaximum("nDiMuons"))];
      _dimuonList  = new IDiMuonList();

      _pDST->SetBranchAddress("dimuon",  &_dimuonList);
    }
  else
    {
      _dimuon     = NULL;
      _dimuonList = NULL;
    }


  // if bbc object exists then add branch
  if (_pDST->GetMaximum("nBbcs") > 0)
    {
      _bbc     = new IBbc();
      _bbcList = new IBbcList();

      _pDST->SetBranchAddress("bbc", &_bbcList);
    }
  else
    {
      _bbc     = NULL;
      _bbcList = NULL;
    }


  // if zdc object exists then add branch
  if (_pDST->GetMaximum("nZdcSmds") > 0)
    {
      _zdcsmd     = new IZdcSmd();
      _zdcsmdList = new IZdcSmdList();

      _pDST->SetBranchAddress("zdcsmd", &_zdcsmdList);
    }
  else
    {
      _zdcsmd     = NULL;
      _zdcsmdList = NULL;
    }
}
//-----------------------------------------------------------------------------


IPdst::~IPdst()
{
  delete _pDST;

  delete _header;

  // all the list objects are most probably deletet by root while delete _pDST
  if (_photon)     delete [] _photon;
  if (_track)      delete [] _track;
  if (_elect)      delete [] _elect;
  if (_muon)       delete [] _muon;
  if (_dimuon)     delete [] _dimuon;
  if (_bbc)        delete _bbc;
  if (_zdcsmd)     delete _zdcsmd;
}
//-----------------------------------------------------------------------------


int IPdst::GetEntries()
{
  return _pDST->GetEntries();
}
//-----------------------------------------------------------------------------


void IPdst::GetEntry(int i)
{
  //_zdcsmd = NULL;

  _pDST->GetEntry(i);

  _header = _headerList->GetHeader();

  if (_bbcList)    _bbc    = _bbcList->GetBbc();
  if (_zdcsmdList)
    {
      _nNeutrons = _zdcsmdList->GetNZdcSmd();
      _zdcsmd    = _zdcsmdList->GetZdcSmd();
    }

  if (_photonList)
    {
      _nPhotons = _photonList->GetNPhotons();

      for (unsigned int i=0; i<_nPhotons; i++)
	_photon[i] = _photonList->GetPhoton(i);
    }


  if (_trackList)
    {
      _nTracks  = _trackList->GetNTracks();

      for (unsigned int i=0; i<_nTracks; i++)
	_track[i] = _trackList->GetTrack(i);
    }

  if (_electList)
    {
      _nElectrons  = _electList->GetNElectrons();

      for (unsigned int i=0; i<_nElectrons; i++)
	_elect[i] = _electList->GetElectron(i);
    }

  if (_muonList)
    {
      _nMuons  = _muonList->GetNMuons();

      for (unsigned int i=0; i<_nMuons; i++)
	_muon[i] = _muonList->GetMuon(i);
    }

  if (_dimuonList)
    {
      _nDiMuons  = _dimuonList->GetNDiMuons();

      for (unsigned int i=0; i<_nDiMuons; i++)
	_dimuon[i] = _dimuonList->GetDiMuon(i);
    }
}
//-----------------------------------------------------------------------------


void IPdst::SetBranchStatus(const char* brname, const int status)
{
  _pDST->SetBranchStatus(brname, status);
}
//-----------------------------------------------------------------------------


bool IPdst::CompareTowerID(short towerID)
{
  bool trackExists = false;

  for (unsigned int i = 0; i < _nTracks; i++)
    if (towerID == _track[i]->GetTowerID())
      {
	trackExists = true;
	break;
      }

  return trackExists;
}
//-----------------------------------------------------------------------------


int IPdst::GetRunNumberFromPDST()
{
  int runNum = -9999;
  for (long evt=0; evt<_pDST->GetEntries(); evt++)
    {
      _pDST->GetEntry(evt);

      if (runNum != _headerList->GetHeader()->GetRunID())
	{
	  runNum = _headerList->GetHeader()->GetRunID();
	  break;
	}
    }
  return runNum;
}
//-----------------------------------------------------------------------------


double IPdst::DiPhotonMass(short i, short j)
{
  TLorentzVector pi0(0., 0., 0., 0.);

  pi0.SetPxPyPzE(_photon[i]->GetPx() + _photon[j]->GetPx(),
		 _photon[i]->GetPy() + _photon[j]->GetPy(),
		 _photon[i]->GetPz() + _photon[j]->GetPz(),
		 _photon[i]->GetEn() + _photon[j]->GetEn());

  return pi0.M();
}
//------------------------------------------------------------------------------


double IPdst::DiHadronMass(short i, short j)
{
  TLorentzVector dihad(0., 0., 0., 0.);

  dihad.SetPxPyPzE(_track[i]->GetPx() + _track[j]->GetPx(),
		   _track[i]->GetPy() + _track[j]->GetPy(),
		   _track[i]->GetPz() + _track[j]->GetPz(),
		   _track[i]->GetP()  + _track[j]->GetP());

  return dihad.M();
}
//------------------------------------------------------------------------------


double IPdst::DiMuonMass(short i, short j)
{
  TLorentzVector dimu(0., 0., 0., 0.);

  dimu.SetPxPyPzE(_muon[i]->GetPx(IMuon::VTX) + _muon[j]->GetPx(IMuon::VTX),
		  _muon[i]->GetPy(IMuon::VTX) + _muon[j]->GetPy(IMuon::VTX),
		  _muon[i]->GetPz(IMuon::VTX) + _muon[j]->GetPz(IMuon::VTX),
		  _muon[i]->GetP (IMuon::VTX) + _muon[j]->GetP (IMuon::VTX));

  return dimu.M();
}
//------------------------------------------------------------------------------


double IPdst::ECone(unsigned short i, double CAng)
{
  double phi = _photon[i]->GetPhi();
  if (phi < -PIby2) phi += TwoPI; // make it consistant with PHENIX convension

  TVector3 vPhoton(0., 0., 0.);
  vPhoton.SetPtThetaPhi(_photon[i]->GetPt(),
			_photon[i]->GetTheta(), phi);

  TVector3 vParticle(0., 0., 0.);

  double econe = 0.;
  double angle = 0.;

  for (unsigned int j=0; j<_nPhotons; j++)
    {
      if (i == j) continue;

      phi = _photon[j]->GetPhi();
      if (phi < -PIby2) phi += TwoPI;

      vParticle.SetPtThetaPhi(_photon[j]->GetPt(),
			      _photon[j]->GetTheta(), phi);

      angle = vPhoton.Angle(vParticle);

      if (angle < CAng)	econe += _photon[j]->GetEnergy();
    }

  for (unsigned int j=0; j<_nTracks; j++)
    {
      vParticle.SetPtThetaPhi(_track[j]->GetPt(),
			      _track[j]->GetTheta0(),
			      _track[j]->GetPhi0());

      angle = vPhoton.Angle(vParticle);

      if (angle < CAng)	econe += _track[j]->GetPt();
    }

  return econe;
}
//------------------------------------------------------------------------------
