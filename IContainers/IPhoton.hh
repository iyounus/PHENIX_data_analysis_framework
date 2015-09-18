//-----------------------------------------------------------
// IPhoton.hh
//	 Created  April 11 2009
//       Imran Younus
//-----------------------------------------------------------

#ifndef IPHOTON_H
#define IPHOTON_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <math.h>



class IPhoton : public TObject
{
private:

  short TowerID;
  short NTowers;           // Multiplicity

  float Energy;
  float X;
  float Y;
  float Z;

  float Prob;
  float Tof;

  short Ert;
  /*
    trigger bit in ErtOut are 
    0:4x4a, 1:4x4b, 2:4x4c, 3:2x2, 4:RICH

    So, here I'm calculating Ert as a 5 digit decimel number:
    RICH*10^4 + 2x2*10^3 + 4x4c*10^2 + 4x4b*10 + 4x4a

    This is how you can get which trigger was fired:
    for example, if you need to check 4x4c:

    (Ert%pow(10, 3) - Ert%pow(10, 2))/Ert%pow(10, 2)
  */


public:
  IPhoton(){;}
  ~IPhoton(){;}


  short  GetNTowers() const {return NTowers;}
  short  GetTowerID() const {return TowerID;}

  // In EMCal convension, West arm is 0 and East is 1.
  // So if TowerID < 4*36*72, then arm is 0
  short  GetArm()     const {return TowerID < 10368 ? 0 : 1;} 
  bool   IsPbSc()     const {return TowerID < 15552 ? true : false;}
  bool   IsPbGl()     const {return TowerID > 15551 ? true : false;}
  short  GetSector();
  void   GetTowerPos(short &sector, short &ypos, short &zpos); 


  float  GetEn()      const {return Energy;}
  float  GetEnergy()  const {return Energy;}
  float  GetX()       const {return X;}
  float  GetY()       const {return Y;}
  float  GetZ()       const {return Z;}

  float  GetProb()    const {return Prob;}
  float  GetTof()     const {return Tof;}

  double GetPhi()     const {return atan2(Y,X);}
  double GetTheta()   const {return atan2(sqrt(X*X + Y*Y), Z);}

  double GetPt()      const {return Energy * sin(GetTheta());}
  double GetPx()      const {return Energy * sin(GetTheta()) * cos(GetPhi());}
  double GetPy()      const {return Energy * sin(GetTheta()) * sin(GetPhi());}
  double GetPz()      const {return Energy * cos(GetTheta());}

  short  GetErt()     const {return Ert;}
  bool   FiredErt()   const {return bool(Ert);}
  bool   FiredRICH()  const {return bool((Ert - Ert%10000)/10000.);}
  bool   Fired2x2()   const {return bool((Ert%10000 - Ert%1000)/1000.);}
  bool   Fired4x4c()  const {return bool((Ert%1000 - Ert%100)/100.);}
  bool   Fired4x4b()  const {return bool((Ert%100 - Ert%10)/10.);}
  bool   Fired4x4a()  const {return bool(Ert%10);}


  void SetTowerID(short id) {TowerID = id;}
  void SetNTowers(short nn) {NTowers = nn;}

  void SetEnergy(float ener) {Energy = ener;}

  void SetX(float x)   {X = x;}
  void SetY(float y)   {Y = y;}
  void SetZ(float z)   {Z = z;}

  void SetProb(float prob) {Prob = prob;}
  void SetTof(float tof)   {Tof = tof;}
  void SetErt(short ert)   {Ert = ert;}

  void Copy(IPhoton *photon);

protected:
  ClassDef(IPhoton,1)
};
//=============================================================================

//inline functions
inline
short IPhoton::GetSector()
{
  if (GetArm() == 0) return int(TowerID/2592); // West arm
  if (GetArm() == 1 && IsPbSc()) return int(TowerID/2592)-2;
  if (IsPbGl()) return int((TowerID-15552)/4608);
  return -9;
}
//-----------------------------------------------------------------------------


inline
void IPhoton::GetTowerPos(short &sect, short &ypos, short &zpos)
{
  short arm = GetArm();
  sect = GetSector();
  int tw = -9;

  if (arm == 0)
    for (int y=0; y<36; y++)
      for (int z=0; z<72; z++)
	{
	  tw = z + 72*y + sect*2592;
	  if (tw == TowerID)
	    {
	      ypos = y;
	      zpos = z;
	      break;
	    }
	}

  if (arm == 1 && sect > 1)
    for (int y=0; y<36; y++)
      for (int z=0; z<72; z++)
	{
	  tw = z + 72*y + (sect+2)*2592;
	  if (tw == TowerID)
	    {
	      ypos = y;
	      zpos = z;
	      break;
	    }
	}

  if (arm == 1 && sect < 2)
    for (int y=0; y<48; y++)
      for (int z=0; z<96; z++)
	{
	  tw = z + 96*y + sect*4608 + 15552;
	  if (tw == TowerID)
	    {
	      ypos = y;
	      zpos = z;
	      break;
	    }
	}
}
//-----------------------------------------------------------------------------


inline
void IPhoton::Copy(IPhoton *photon)
{
  TowerID  = photon->TowerID;
  NTowers  = photon->NTowers;

  Energy   = photon->Energy;
  X        = photon->X;
  Y        = photon->Y;
  Z        = photon->Z;

  Prob     = photon->Prob;
  Tof      = photon->Tof;

  Ert      = photon->Ert;
}
//-----------------------------------------------------------------------------
#endif
