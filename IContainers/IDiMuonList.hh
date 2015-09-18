//-----------------------------------------------------------
// IDiMuonList.hh
//       Created Feb 23 2011
//       Imran Younus
//-----------------------------------------------------------

#ifndef IDIMUONLIST_HH
#define IDIMUONLIST_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TClonesArray;
class IDiMuon;

class IDiMuonList : public TObject
{
public:
  IDiMuonList();
  virtual ~IDiMuonList();

  void Reset();

  unsigned short GetNDiMuons() const { return nDiMuons; }

  IDiMuon* AddDiMuon();
  IDiMuon* GetDiMuon(const unsigned int idimu);

protected:
  unsigned short nDiMuons;
  TClonesArray *dimuonList;

private:
  ClassDef(IDiMuonList,1);
};
#endif
