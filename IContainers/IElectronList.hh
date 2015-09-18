//-----------------------------------------------------------
// IElectronList.hh
//       Created Feb 09, 2010
//       Imran Younus
//-----------------------------------------------------------

#ifndef IELECTRONLIST_HH
#define IELECTRONLIST_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TClonesArray;
class IElectron;

class IElectronList : public TObject
{
public:
  IElectronList();
  virtual ~IElectronList();

  void Reset();

  unsigned short GetNElectrons() const { return nElectrons; }

  IElectron* AddElectron();
  IElectron* GetElectron(const unsigned int ielectron);

protected:
  unsigned short nElectrons;
  TClonesArray *electronList;

private:
  ClassDef(IElectronList,1);
};
#endif
