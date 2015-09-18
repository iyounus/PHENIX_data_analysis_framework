#include "IMuonList.hh"

#include "TClonesArray.h"
#include "IMuon.hh"


ClassImp(IMuonList)

#define NMUONS 200

IMuonList::IMuonList()
{
  nMuons = 0;
  muonList = new TClonesArray("IMuon", NMUONS);
}
//-----------------------------------------------------------------------------


IMuonList::~IMuonList()
{
  muonList->Clear();
  delete muonList;
}
//-----------------------------------------------------------------------------


void IMuonList::Reset()
{
  muonList->Clear();
  nMuons = 0;
}
//-----------------------------------------------------------------------------


IMuon* IMuonList::AddMuon()
{
  TClonesArray &plist = *muonList;
  IMuon *muon = new(plist[nMuons++]) IMuon();
  return muon;
}
//-----------------------------------------------------------------------------


IMuon* IMuonList::GetMuon(const unsigned int itrk)
{
  return (IMuon*)muonList->UncheckedAt(itrk);
}
//-----------------------------------------------------------------------------
