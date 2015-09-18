//-----------------------------------------------------------
// IPhotonList.hh
//       Created April 11 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IPHOTONLIST_HH
#define IPHOTONLIST_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TClonesArray;
class IPhoton;

class IPhotonList : public TObject
{
public:
  IPhotonList();
  virtual ~IPhotonList();

  void Reset();

  unsigned short GetNPhotons() const { return nPhotons; }

  IPhoton* AddPhoton();
  IPhoton* GetPhoton(const unsigned int iphoton);

protected:
  unsigned short nPhotons;
  TClonesArray *photonList;

private:
  ClassDef(IPhotonList,1);
};
#endif
