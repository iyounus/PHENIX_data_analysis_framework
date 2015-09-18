//-----------------------------------------------------------
// ITrack.hh
//	 Created  April 11 2009
//       Imran Younus
//-----------------------------------------------------------
#ifndef ITRACK_H
#define ITRACK_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <math.h>

class ITrack : public TObject
{
private:
  short Charge;       // PHCentralTrack::get_charge
  short Quality;      // track quality (>3)
  short TowerID;      // EMCal tower ID corresponding to track

  float P;            // PHCentralTrack::get_mom
  float Phi0;         // PHCentralTrack::get_phi0
  float Theta0;       // PHCentralTrack::get_the0
  float PhiDC;        // PHCentralTrack::get_phi
  float ZedDC;        // PHCentralTrack::get_zed
  float Alpha;        // PHCentralTrack::get_alpha

  float Pc3sdPhi;     // PHCentralTrack::get_pc3sdphi
  float Pc3sdZ;       // PHCentralTrack::get_pc3sdz

  float EmcsdPhi;     // PHCentralTrack::get_emcsdphi
  float EmcsdZ;       // PHCentralTrack::get_emcsdz

  float Ecore;        // PHCentralTrack::get_ecore
  float EmcE;         // PHCentralTrack::get_emce
  float Prob;         // PHCentralTrack::get_prob

  short Ert;
  /*
    trigger bit in ErtOut are 
    0:4x4a, 1:4x4b, 2:4x4c, 3:2x2, 4:RICH

    So, here I'm calculating Ert as a 5 digit decimel number:

    RICH*10^4 + 2x2*10^3 + 4x4c*10^2 + 4x4b*10 + 4x4a

    This is how you can get which trigger was fired:
    for example, if you need to chech 4x4c:

    (Ert%pow(10, 3) - Ert%pow(10, 2))/Ert%pow(10, 2)
  */


public:
  ITrack(){;}
  ~ITrack(){;}


  short  GetCharge()  const {return Charge;}
  short  GetQuality() const {return Quality;}
  short  GetTowerID() const {return TowerID;}

  short  GetErt()     const {return Ert;}
  bool   FiredErt()   const {return bool(Ert);}
  bool   FiredRICH()  const {return bool((Ert - Ert%10000)/10000.);}
  bool   Fired2x2()   const {return bool((Ert%10000 - Ert%1000)/1000.);}
  bool   Fired4x4c()  const {return bool((Ert%1000 - Ert%100)/100.);}
  bool   Fired4x4b()  const {return bool((Ert%100 - Ert%10)/10.);}
  bool   Fired4x4a()  const {return bool(Ert%10);}

  double GetPt()      const {return P * sin(Theta0);}
  double GetPx()      const {return P * sin(Theta0) * cos(Phi0);}
  double GetPy()      const {return P * sin(Theta0) * sin(Phi0);}
  double GetPz()      const {return P * cos(Theta0);}


  float  GetP()        const {return P;}
  float  GetMomentum() const {return P;}
  float  GetPhi0()     const {return Phi0;}
  float  GetTheta0()   const {return Theta0;}
  float  GetPhiDC()    const {return PhiDC;}
  float  GetZedDC()    const {return ZedDC;}
  float  GetAlpha()    const {return Alpha;}


  float  GetPc3sdPhi() const {return Pc3sdPhi;}
  float  GetPc3sdZ()   const {return Pc3sdZ;}

  double GetPc3MatchingSig()
  {return sqrt(Pc3sdPhi*Pc3sdPhi + Pc3sdZ*Pc3sdZ);}

  float GetEmcsdPhi() const {return EmcsdPhi;}
  float GetEmcsdZ()   const {return EmcsdZ;}

  double GetEmcMatchingSig()
  {return sqrt(EmcsdPhi*EmcsdPhi + EmcsdZ*EmcsdZ);}

  float GetEcore() const {return Ecore;}
  float GetEmcE()  const {return EmcE;}
  float GetProb()  const {return Prob;}

  float GetEPratio() const {return EmcE/P;}
  float GetEoverP()  const {return EmcE/P;}

  void SetCharge(short q)   {Charge = q;}
  void SetQuality(short q)  {Quality = q;}
  void SetTowerID(short id) {TowerID = id;}
  void SetErt(short ert)    {Ert = ert;}


  void SetMomentum(float p) {P = p;}
  void SetPhi0(float ph)    {Phi0 = ph;}
  void SetTheta0(float th)  {Theta0 = th;}
  void SetPhiDC(float ph)   {PhiDC = ph;}
  void SetZedDC(float zed)  {ZedDC = zed;}
  void SetAlpha(float alp)  {Alpha = alp;}


  void SetPc3sdPhi(float sdPhi) {Pc3sdPhi = sdPhi;}
  void SetPc3sdZ  (float sdZ)   {Pc3sdZ = sdZ;}

  void SetEmcsdPhi(float sdPhi) {EmcsdPhi = sdPhi;}
  void SetEmcsdZ  (float sdZ)   {EmcsdZ = sdZ;}

  void SetEcore(float e)   {Ecore = e;}
  void SetEmcE(float emce) {EmcE = emce;}
  void SetProb(float prob) {Prob = prob;}


  void Copy(ITrack *track);


protected:
  ClassDef(ITrack,1)
};
//=============================================================================


// inline functions

inline
void ITrack::Copy(ITrack *track)
{
  Charge   = track->Charge;
  Quality  = track->Quality;
  TowerID  = track->TowerID;
  Ert      = track->Ert;

  P        = track->P;
  Phi0     = track->Phi0;
  Theta0   = track->Theta0;
  PhiDC    = track->PhiDC;
  ZedDC    = track->ZedDC;
  Alpha    = track->Alpha;

  Pc3sdPhi = track->Pc3sdPhi;
  Pc3sdZ   = track->Pc3sdZ;

  EmcsdPhi = track->EmcsdPhi;
  EmcsdZ   = track->EmcsdZ;

  Ecore    = track->Ecore;
  EmcE     = track->EmcE;
  Prob     = track->Prob;
}
//-----------------------------------------------------------------------------
#endif
