#include "IPhotonList.hh"

#include "TClonesArray.h"
#include "IPhoton.hh"


ClassImp(IPhotonList)

#define NPHOTONS 200

IPhotonList::IPhotonList()
{
  nPhotons = 0;
  photonList = new TClonesArray("IPhoton", NPHOTONS);
}
//-----------------------------------------------------------------------------


IPhotonList::~IPhotonList()
{
  photonList->Clear();
  delete photonList;
}
//-----------------------------------------------------------------------------


void IPhotonList::Reset()
{
  photonList->Clear();
  nPhotons = 0;
}
//-----------------------------------------------------------------------------


IPhoton* IPhotonList::AddPhoton()
{
  TClonesArray &plist = *photonList;
  IPhoton *photon = new(plist[nPhotons++]) IPhoton();
  return photon;
}
//-----------------------------------------------------------------------------


IPhoton* IPhotonList::GetPhoton(const unsigned int itrk)
{
  return (IPhoton*)photonList->UncheckedAt(itrk);
}
//-----------------------------------------------------------------------------
