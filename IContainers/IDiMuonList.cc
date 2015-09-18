#include "IDiMuonList.hh"

#include "TClonesArray.h"
#include "IDiMuon.hh"


ClassImp(IDiMuonList)

#define NDIMUONS 200

IDiMuonList::IDiMuonList()
{
  nDiMuons = 0;
  dimuonList = new TClonesArray("IDiMuon", NDIMUONS);
}
//-----------------------------------------------------------------------------


IDiMuonList::~IDiMuonList()
{
  dimuonList->Clear();
  delete dimuonList;
}
//-----------------------------------------------------------------------------


void IDiMuonList::Reset()
{
  dimuonList->Clear();
  nDiMuons = 0;
}
//-----------------------------------------------------------------------------


IDiMuon* IDiMuonList::AddDiMuon()
{
  TClonesArray &plist = *dimuonList;
  IDiMuon *dimuon = new(plist[nDiMuons++]) IDiMuon();
  return dimuon;
}
//-----------------------------------------------------------------------------


IDiMuon* IDiMuonList::GetDiMuon(const unsigned int itrk)
{
  return (IDiMuon*)dimuonList->UncheckedAt(itrk);
}
//-----------------------------------------------------------------------------
