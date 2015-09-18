#include "IEvent.hh"

//#include <iostream>
//using namespace std;

#include "IHeader.hh"
#include "IBbc.hh"
#include "IZdcSmd.hh"
#include "IPhoton.hh"
#include "ITrack.hh"
#include "IElectron.hh"
#include "IMuon.hh"
#include "IDiMuon.hh"


#include "IConsts.hh"

#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;

IEvent::IEvent()
{
  _header = NULL;
  _bbc    = NULL;   // bbc and zdcsmd re added only when used;
  _zdcsmd = NULL;

  // these vectors may or maynot be used by the user
  _photon = new vector<IPhoton *>;
  _track  = new vector<ITrack *>;
  _elect  = new vector<IElectron *>;
  _muon   = new vector<IMuon *>;
  _dimu   = new vector<IDiMuon *>;
}
//------------------------------------------------------------------------------


IEvent::~IEvent()
{
  if (_header ) delete _header;
  if (_bbc)     delete _bbc;
  if (_zdcsmd)  delete _zdcsmd;


  for (unsigned int i=0; i<_photon->size(); i++)
    delete _photon->at(i);

  _photon->clear();
  delete _photon;

  for (unsigned int i=0; i<_track->size(); i++)
    delete _track->at(i);

  _track->clear();
  delete _track;

  for (unsigned int i=0; i<_elect->size(); i++)
    delete _elect->at(i);

  _elect->clear();
  delete _elect;

  for (unsigned int i=0; i<_muon->size(); i++)
    delete _muon->at(i);

  _muon->clear();
  delete _muon;

  for (unsigned int i=0; i<_dimu->size(); i++)
    delete _dimu->at(i);

  _dimu->clear();
  delete _dimu;
}
//-----------------------------------------------------------------------------


void IEvent::AddHeader(IHeader  *hdr)
{
  _header = new IHeader();
  _header->Copy(hdr);
}
//-----------------------------------------------------------------------------


void IEvent::AddBbc(IBbc *bbc)
{
  _bbc = new IBbc();
  _bbc->Copy(bbc);
}
//-----------------------------------------------------------------------------


void IEvent::AddZdcSmd(IZdcSmd *zdc)
{
  _zdcsmd = new IZdcSmd();
  _zdcsmd->Copy(zdc);
}
//-----------------------------------------------------------------------------


void IEvent::AddPhoton(IPhoton *pho)
{
  IPhoton *photon = new IPhoton();
  photon->Copy(pho);
  _photon->push_back(photon);
}
//-----------------------------------------------------------------------------


void IEvent::AddTrack (ITrack *trk)
{
  ITrack *track = new ITrack();
  track->Copy(trk);
  _track->push_back(track);
}
//-----------------------------------------------------------------------------


void IEvent::AddElectron(IElectron *ele)
{
  IElectron *elect = new IElectron();
  elect->Copy(ele);
  _elect->push_back(elect);
}
//-----------------------------------------------------------------------------


void IEvent::AddMuon(IMuon *mu)
{
  IMuon *muon = new IMuon();
  muon->Copy(mu);
  _muon->push_back(muon);
}
//-----------------------------------------------------------------------------


void IEvent::AddDiMuon(IDiMuon *dimu)
{
  IDiMuon *dimuon = new IDiMuon();
  dimuon->Copy(dimu);
  _dimu->push_back(dimuon);
}
//-----------------------------------------------------------------------------


double IEvent::DiPhotonAngle(short i, short j)
{
  IPhoton *photon1 = _photon->at(i);
  IPhoton *photon2 = _photon->at(j);

  TVector3 v1(photon1->GetPx(), photon1->GetPy(), photon1->GetPz());
  TVector3 v2(photon2->GetPx(), photon2->GetPy(), photon2->GetPz());

  return v1.Angle(v2);
}
//-----------------------------------------------------------------------------


double IEvent::DiPhotonMass(short i, short j)
{
  TLorentzVector pi0(0., 0., 0., 0.);

  IPhoton *photon1 = _photon->at(i);
  IPhoton *photon2 = _photon->at(j);

  pi0.SetPxPyPzE(photon1->GetPx() + photon2->GetPx(),
		 photon1->GetPy() + photon2->GetPy(),
		 photon1->GetPz() + photon2->GetPz(),
		 photon1->GetEn() + photon2->GetEn());

  return pi0.M();
}
//------------------------------------------------------------------------------


double IEvent::DiHadronMass(short i, short j)
{
  TLorentzVector dihad(0., 0., 0., 0.);

  ITrack *track1 = _track->at(i);
  ITrack *track2 = _track->at(j);

  dihad.SetPxPyPzE(track1->GetPx() + track2->GetPx(),
		   track1->GetPy() + track2->GetPy(),
		   track1->GetPz() + track2->GetPz(),
		   track1->GetP()  + track2->GetP());

  return dihad.M();
}
//------------------------------------------------------------------------------


double IEvent::DiMuonMass(short i, short j)
{
  TLorentzVector dimu(0., 0., 0., 0.);

  IMuon *muon1 = _muon->at(i);
  IMuon *muon2 = _muon->at(j);

  dimu.SetPxPyPzE(muon1->GetPx(IMuon::VTX) + muon2->GetPx(IMuon::VTX),
		  muon1->GetPy(IMuon::VTX) + muon2->GetPy(IMuon::VTX),
		  muon1->GetPz(IMuon::VTX) + muon2->GetPz(IMuon::VTX),
		  muon1->GetP (IMuon::VTX) + muon2->GetP (IMuon::VTX));

  return dimu.M();
}
//------------------------------------------------------------------------------


double IEvent::ECone(unsigned short i, double CAng)
{
  double phi = _photon->at(i)->GetPhi();
  if (phi < -PIby2) phi += TwoPI; // make it consistant with PHENIX convension

  TVector3 vPhoton(0., 0., 0.);
  vPhoton.SetPtThetaPhi(_photon->at(i)->GetPt(),
			_photon->at(i)->GetTheta(), phi);

  TVector3 vParticle(0., 0., 0.);

  double econe = 0.;
  double angle = 0.;

  for (unsigned int j=0; j<_photon->size(); j++)
    {
      if (i == j) continue;

      phi = _photon->at(j)->GetPhi();
      if (phi < -PIby2) phi += TwoPI;

      vParticle.SetPtThetaPhi(_photon->at(j)->GetPt(),
			      _photon->at(j)->GetTheta(), phi);

      angle = vPhoton.Angle(vParticle);

      if (angle < CAng)	econe += _photon->at(j)->GetEnergy();
    }

  for (unsigned int j=0; j<_track->size(); j++)
    {
      vParticle.SetPtThetaPhi(_track->at(j)->GetPt(),
			      _track->at(j)->GetTheta0(),
			      _track->at(j)->GetPhi0());

      angle = vPhoton.Angle(vParticle);

      if (angle < CAng)	econe += _track->at(j)->GetPt();
    }

  return econe;
}
//------------------------------------------------------------------------------
