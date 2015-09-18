//-----------------------------------------------------------
// IZdcSmdList.hh
//       Created April 11 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IZDCSMDLIST_HH
#define IZDCSMDLIST_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TClonesArray;
class IZdcSmd;

class IZdcSmdList : public TObject
{
public:
  IZdcSmdList();
  virtual ~IZdcSmdList();

  void Reset();

  unsigned short GetNZdcSmd() const { return nZdcSmds; }

  IZdcSmd* AddZdcSmd();
  IZdcSmd* GetZdcSmd();

protected:
  unsigned short nZdcSmds;
  TClonesArray *zdcsmdList;

private:
  ClassDef(IZdcSmdList,1);
};
#endif
