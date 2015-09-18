//-----------------------------------------------------------
// IHeaderList.hh
//       Created Aug 06 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IBBCLIST_HH
#define IBBCLIST_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TClonesArray;
class IBbc;

class IBbcList : public TObject
{
public:
  IBbcList();
  virtual ~IBbcList();

  void Reset();

  unsigned short GetNBbc() const { return nBbcs; }

  IBbc* AddBbc();
  IBbc* GetBbc();

protected:
  unsigned short nBbcs;
  TClonesArray *bbcList;

private:
  ClassDef(IBbcList,1);
};
#endif
