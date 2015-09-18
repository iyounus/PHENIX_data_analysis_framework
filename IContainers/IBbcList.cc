#include "IBbcList.hh"

#include "TClonesArray.h"
#include "IBbc.hh"


ClassImp(IBbcList)

#define NBBCS 1 // one bbc object per event

IBbcList::IBbcList()
{
  nBbcs = 0;
  bbcList = new TClonesArray("IBbc", NBBCS);
}
//-----------------------------------------------------------------------------


IBbcList::~IBbcList()
{
  bbcList->Clear();
  delete bbcList;
}
//-----------------------------------------------------------------------------


void IBbcList::Reset()
{
  bbcList->Clear();
  nBbcs = 0;
}
//-----------------------------------------------------------------------------


IBbc* IBbcList::AddBbc()
{
  // there is only one bbc per event
  if (nBbcs>0) return NULL;
  TClonesArray &hlist = *bbcList;
  IBbc *bbc = new(hlist[nBbcs++]) IBbc();
  return bbc;
}
//-----------------------------------------------------------------------------


IBbc* IBbcList::GetBbc()
{
  return (IBbc*)bbcList->UncheckedAt(0);
}
//-----------------------------------------------------------------------------
