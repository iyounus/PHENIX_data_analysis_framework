//-----------------------------------------------------------
// IZdcSmd.hh
//       Created April 11 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IZDCSMD_HH
#define IZDCSMD_HH

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class IZdcSmd : public TObject
{
private:

  float ZdcN_E;
  float ZdcS_E;

  float SmdN_E;
  float SmdS_E;
  float SmdN_X;
  float SmdS_X;
  float SmdN_Y;
  float SmdS_Y;

  float SmdN_chargeY[8];
  float SmdN_chargeX[7];

  float SmdS_chargeY[8];
  float SmdS_chargeX[7];

  unsigned short ScintN_ADC;
  unsigned short ScintS_ADC;

  unsigned short ScintN_TDC;
  unsigned short ScintS_TDC;

 public:

  IZdcSmd(){;}              // constructor
  virtual ~IZdcSmd(){;}     // destructor

  float GetZdcN_E() const {return ZdcN_E;}
  float GetZdcS_E() const {return ZdcS_E;}

  float GetSmdN_E() const {return SmdN_E;}
  float GetSmdS_E() const {return SmdS_E;}
  float GetSmdN_X() const {return SmdN_X;}
  float GetSmdS_X() const {return SmdS_X;}
  float GetSmdN_Y() const {return SmdN_Y;}
  float GetSmdS_Y() const {return SmdS_Y;}

  float GetSmdN_chargeY(int ipmt) {return  SmdN_chargeY[ipmt];}
  float GetSmdN_chargeX(int ipmt) {return  SmdN_chargeX[ipmt];}
  float GetSmdS_chargeY(int ipmt) {return  SmdS_chargeY[ipmt];}
  float GetSmdS_chargeX(int ipmt) {return  SmdS_chargeX[ipmt];}

  unsigned short GetScintN_ADC() const {return ScintN_ADC;}
  unsigned short GetScintS_ADC() const {return ScintS_ADC;}
  unsigned short GetScintN_TDC() const {return ScintN_TDC;}
  unsigned short GetScintS_TDC() const {return ScintS_TDC;}

  void SetZdcN_E(float ne) {ZdcN_E = ne;}
  void SetZdcS_E(float se) {ZdcS_E = se;}

  void SetSmdN_E(float ne) {SmdN_E = ne;}
  void SetSmdS_E(float se) {SmdS_E = se;}
  void SetSmdN_X(float nx) {SmdN_X = nx;}
  void SetSmdS_X(float sx) {SmdS_X = sx;}
  void SetSmdN_Y(float ny) {SmdN_Y = ny;}
  void SetSmdS_Y(float sy) {SmdS_Y = sy;}

  void SetSmdN_chargeY(int i, float ch) {SmdN_chargeY[i] = ch;}
  void SetSmdN_chargeX(int i, float ch) {SmdN_chargeX[i] = ch;}
  void SetSmdS_chargeY(int i, float ch) {SmdS_chargeY[i] = ch;}
  void SetSmdS_chargeX(int i, float ch) {SmdS_chargeX[i] = ch;}

  void SetSmdN_chargeY(float *ch)
  {for (int i=0; i<8; i++) SmdN_chargeY[i] = ch[i];}

  void SetSmdN_chargeX(float *ch)
  {for (int i=0; i<7; i++) SmdN_chargeX[i] = ch[i];}

  void SetSmdS_chargeY(float *ch)
  {for (int i=0; i<8; i++) SmdS_chargeY[i] = ch[i];}

  void SetSmdS_chargeX(float *ch)
  {for (int i=0; i<7; i++) SmdS_chargeX[i] = ch[i];}

  void SetScintN_ADC(unsigned short adc) {ScintN_ADC = adc;}
  void SetScintS_ADC(unsigned short adc) {ScintS_ADC = adc;}
  void SetScintN_TDC(unsigned short tdc) {ScintN_TDC = tdc;}
  void SetScintS_TDC(unsigned short tdc) {ScintS_TDC = tdc;}

  void Copy(IZdcSmd *zdc);


protected:
  ClassDef(IZdcSmd,1)
};
//=============================================================================

//inline functions


inline
void IZdcSmd::Copy(IZdcSmd *zdc)
{
  ZdcN_E = zdc->ZdcN_E;
  ZdcS_E = zdc->ZdcS_E;

  SmdN_E = zdc->SmdN_E;
  SmdS_E = zdc->SmdS_E;
  SmdN_X = zdc->SmdN_X;
  SmdS_X = zdc->SmdS_Y;
  SmdN_Y = zdc->SmdN_Y;
  SmdS_Y = zdc->SmdS_Y;


  for (int i=0; i<8; i++)
    {
      SmdN_chargeY[i] = zdc->SmdN_chargeY[i];
      SmdS_chargeY[i] = zdc->SmdS_chargeY[i];
    }

  for (int i=0; i<7; i++)
    {
      SmdN_chargeX[i] = zdc->SmdN_chargeX[i];
      SmdS_chargeX[i] = zdc->SmdS_chargeX[i];
    }

  ScintN_ADC = zdc->ScintN_ADC;
  ScintS_ADC = zdc->ScintS_ADC;
  ScintN_TDC = zdc->ScintN_TDC;
  ScintS_TDC = zdc->ScintS_TDC;
}
//-----------------------------------------------------------------------------
#endif
