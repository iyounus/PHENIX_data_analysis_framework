#include "IElectronList.hh"

#include "TClonesArray.h"
#include "IElectron.hh"


ClassImp(IElectronList)

#define NELECTRONS 100

IElectronList::IElectronList()
{
  nElectrons = 0;
  electronList = new TClonesArray("IElectron", NELECTRONS);
}
//-----------------------------------------------------------------------------


IElectronList::~IElectronList()
{
  electronList->Clear();
  delete electronList;
}
//-----------------------------------------------------------------------------


void IElectronList::Reset()
{
  electronList->Clear();
  nElectrons = 0;
}
//-----------------------------------------------------------------------------


IElectron* IElectronList::AddElectron()
{
  TClonesArray &elist = *electronList;
  IElectron *electron = new(elist[nElectrons++]) IElectron();
  return electron;
}
//-----------------------------------------------------------------------------


IElectron* IElectronList::GetElectron(const unsigned int iele)
{
  return (IElectron*)electronList->UncheckedAt(iele);
}
//-----------------------------------------------------------------------------
