#include "IHeaderList.hh"

#include "TClonesArray.h"
#include "IHeader.hh"


ClassImp(IHeaderList)

#define NHEADERS 1 // one header per event

IHeaderList::IHeaderList()
{
  nHeaders = 0;
  headerList = new TClonesArray("IHeader", NHEADERS);
}
//-----------------------------------------------------------------------------


IHeaderList::~IHeaderList()
{
  headerList->Clear();
  delete headerList;
}
//-----------------------------------------------------------------------------


void IHeaderList::Reset()
{
  headerList->Clear();
  nHeaders = 0;
}
//-----------------------------------------------------------------------------


IHeader* IHeaderList::AddHeader()
{
  // there is only one header per event
  if (nHeaders>0) return NULL;
  TClonesArray &hlist = *headerList;
  IHeader *header = new(hlist[nHeaders++]) IHeader();
  return header;
}
//-----------------------------------------------------------------------------


IHeader* IHeaderList::GetHeader()
{
  return (IHeader*)headerList->UncheckedAt(0);
}
//-----------------------------------------------------------------------------
