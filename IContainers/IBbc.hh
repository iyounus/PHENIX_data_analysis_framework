//-----------------------------------------------------------
// IBbc.hh
//       Created Aug 06 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IBBC_HH
#define IBBC_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class IBbc : public TObject
{
private:

  float BbcT0;

  short BbcMultN;
  short BbcMultS;

  float BbcChargeN;
  float BbcChargeS;


public:
  IBbc(){;}
  ~IBbc(){;}


  float GetBbcT0()   const {return BbcT0;}

  short GetBbcMultN()   const {return BbcMultN;}
  short GetBbcMultS()   const {return BbcMultS;}

  float GetBbcChargeN() const {return BbcChargeN;}
  float GetBbcChargeS() const {return BbcChargeS;}


  void SetBbcT0(float val) {BbcT0 = val;}

  void SetBbcMultN(short multiplicity) {BbcMultN = multiplicity;}
  void SetBbcMultS(short multiplicity) {BbcMultS = multiplicity;}

  void SetBbcChargeN(float charge) {BbcChargeN = charge;}
  void SetBbcChargeS(float charge) {BbcChargeS = charge;}

  void Copy(IBbc *bbc);

protected:
  ClassDef(IBbc,1)
};
//=============================================================================

//inline functions


inline
void IBbc::Copy(IBbc *bbc)
{
  BbcT0 = bbc->GetBbcT0();

  BbcMultN = bbc->GetBbcMultN();
  BbcMultS = bbc->GetBbcMultS();

  BbcChargeN = bbc->GetBbcChargeN();
  BbcChargeS = bbc->GetBbcChargeS();
}
//-----------------------------------------------------------------------------
#endif
