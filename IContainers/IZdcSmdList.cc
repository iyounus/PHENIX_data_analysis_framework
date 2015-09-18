#include "IZdcSmdList.hh"

#include "TClonesArray.h"
#include "IZdcSmd.hh"


ClassImp(IZdcSmdList)

#define NZDCSMD 1 // one zdcsmd per event

IZdcSmdList::IZdcSmdList()
{
  nZdcSmds = 0;
  zdcsmdList = new TClonesArray("IZdcSmd", NZDCSMD);
}
//-----------------------------------------------------------------------------


IZdcSmdList::~IZdcSmdList()
{
  zdcsmdList->Clear();
  delete zdcsmdList;
}
//-----------------------------------------------------------------------------


void IZdcSmdList::Reset()
{
  zdcsmdList->Clear();
  nZdcSmds = 0;
}
//-----------------------------------------------------------------------------


IZdcSmd* IZdcSmdList::AddZdcSmd()
{
  // there is only one zdcsmd per event
  if (nZdcSmds > 0) return NULL;
  TClonesArray &hlist = *zdcsmdList;
  IZdcSmd *zdcsmd = new(hlist[nZdcSmds++]) IZdcSmd();
  return zdcsmd;
}
//-----------------------------------------------------------------------------


IZdcSmd* IZdcSmdList::GetZdcSmd()
{
  return (IZdcSmd*)zdcsmdList->UncheckedAt(0);
}
//-----------------------------------------------------------------------------
