//-----------------------------------------------------------
// IHeaderList.hh
//       Created April 11 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IHEADERLIST_HH
#define IHEADERLIST_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TClonesArray;
class IHeader;

class IHeaderList : public TObject
{
public:
  IHeaderList();
  virtual ~IHeaderList();

  void Reset();

  unsigned short GetNHeaders() const { return nHeaders; }

  IHeader* AddHeader();
  IHeader* GetHeader();

protected:
  unsigned short nHeaders;
  TClonesArray *headerList;

private:
  ClassDef(IHeaderList,1);
};
#endif
