//-----------------------------------------------------------
// IMuonList.hh
//       Created November 15 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IMUONLIST_HH
#define IMUONLIST_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TClonesArray;
class IMuon;

class IMuonList : public TObject
{
public:
  IMuonList();
  virtual ~IMuonList();

  void Reset();

  unsigned short GetNMuons() const { return nMuons; }

  IMuon* AddMuon();
  IMuon* GetMuon(const unsigned int imuon);

protected:
  unsigned short nMuons;
  TClonesArray *muonList;

private:
  ClassDef(IMuonList,1);
};
#endif
