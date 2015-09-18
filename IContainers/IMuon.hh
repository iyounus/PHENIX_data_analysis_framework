//-----------------------------------------------------------
// IMuon.hh
//	 Created  Novemver 15 2009
//       Modified Feb 26, 2011
//       Imran Younus
//-----------------------------------------------------------

#ifndef IMUON_H
#define IMUON_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <math.h>

#define nPOSITIONS 6 // vtx, st1, st2, st3, kalman projection to gap0, gap0

#define MuMass   0.105658367     // from PDG
#define MuMass2  0.0111636905171 // square of muon mass


class IMuon : public TObject
{
private:

  float X[nPOSITIONS]; // vtx, st1, st2, st3, kalman proj, gap0
  float Y[nPOSITIONS];
  float Z[nPOSITIONS];

  float Px[nPOSITIONS-1]; // no momentum at gap0
  float Py[nPOSITIONS-1];
  float Pz[nPOSITIONS-1];

  float Dx;  // dxdz at MuId gap0
  float Dy;  // dydz at MuId gap0

  short Charge;

  unsigned short MuTrHitPat; // hit pattern in MuTr
  unsigned short MuIdHitPat; // hit pattern in MuId

  float MuTrChi2;
  float MuIdChi2;
  float VtxRefitChi2;

public:
  IMuon(){;}
  ~IMuon(){;}

  // vtx, st1, st2, st3, kalman projection to gap0, gap0
  enum  POS {VTX, ST1, ST2, ST3, KAL, GP0};

  float GetX(short pos)   const {return X[pos];}
  float GetY(short pos)   const {return Y[pos];}
  float GetZ(short pos)   const {return Z[pos];}

  float GetPx(short pos)  const {return Px[pos];}
  float GetPy(short pos)  const {return Py[pos];}
  float GetPz(short pos)  const {return Pz[pos];}

  float GetDx()           const {return Dx;};
  float GetDy()           const {return Dy;};
  float GetMuIdSlopeX()   const {return Dx;};
  float GetMuIdSlopeY()   const {return Dy;};

  short GetCharge()       const {return Charge;}

  unsigned short GetMuTrHitPat() const {return MuTrHitPat;}
  unsigned short GetMuIdHitPat() const {return MuIdHitPat;}

  float GetMuTrChi2()     const {return MuTrChi2;}
  float GetMuIdChi2()     const {return MuIdChi2;}

  float GetVtxRefitChi2() const {return VtxRefitChi2;}

  double GetE()           const
  {return sqrt(Px[VTX]*Px[VTX] + Py[VTX]*Py[VTX] + Pz[VTX]*Pz[VTX] + MuMass2);}

  double GetP(short pos)  const
  {return sqrt(Px[pos]*Px[pos] + Py[pos]*Py[pos] + Pz[pos]*Pz[pos]);}

  double GetPt(short pos) const
  {return sqrt(Px[pos]*Px[pos] + Py[pos]*Py[pos]);}

  double GetELoss()       const {return GetP(VTX) - GetP(ST1);}

  double GetEta() const
  {return 0.5*log( (GetE() + Pz[VTX])/(GetE() - Pz[VTX]) );}
  double GetRapidity() const
  {return 0.5*log( (GetE() + Pz[VTX])/(GetE() - Pz[VTX]) );}

  double GetPhi(short pos)   const {return atan2(Py[pos], Px[pos]);}
  double GetTheta(short pos) const {return atan2(GetPt(pos), Pz[pos]);}


  // inline functions
  short  GetNMuTrHits();
  short  GetNMuIdHits();
  double GetDG0();
  double GetDDG0();
  double GetDS3();
  double GetDS3cpt();


  void SetX(short pos, float x)    {X[pos] = x;}
  void SetY(short pos, float y)    {Y[pos] = y;}
  void SetZ(short pos, float z)    {Z[pos] = z;}

  void SetPx(short pos, float px)  {Px[pos] = px;}
  void SetPy(short pos, float py)  {Py[pos] = py;}
  void SetPz(short pos, float pz)  {Pz[pos] = pz;}

  void SetDx(float dx)             {Dx = dx;}
  void SetDy(float dy)             {Dy = dy;}

  void SetCharge(short ch)         {Charge = ch;}

  void SetMuTrHitPat(unsigned short n) {MuTrHitPat = n;}
  void SetMuIdHitPat(unsigned short n) {MuIdHitPat = n;}

  void SetMuTrChi2(float chi)      {MuTrChi2 = chi;}
  void SetMuIdChi2(float chi)      {MuIdChi2 = chi;}

  void SetVtxRefitChi2(float chi2) {VtxRefitChi2 = chi2;}

  void Copy(IMuon *muon);

protected:
  ClassDef(IMuon, 2)
};
//=============================================================================


// inline functions

inline
double IMuon::GetDG0()
{
  double Xproj = 9.e12;
  double Yproj = 9.e12;

  if (X[KAL] == 0 && Y[KAL] == 0)
    {
      Xproj = double(X[ST3] + Px[ST3]* (Z[GP0] - Z[ST3]) /Pz[ST3]);
      Yproj = double(Y[ST3] + Py[ST3]* (Z[GP0] - Z[ST3]) /Pz[ST3]);
    }
  else
    {
      Xproj = double(X[KAL] + Px[ST3]* (Z[GP0] - Z[KAL]) /Pz[ST3]);
      Yproj = double(Y[KAL] + Py[ST3]* (Z[GP0] - Z[KAL]) /Pz[ST3]);
    }


  double dg0 = sqrt( (Xproj-X[GP0])*(Xproj-X[GP0]) +
		     (Yproj-Y[GP0])*(Yproj-Y[GP0]) );
  return dg0;
}
//-----------------------------------------------------------------------------


inline
double IMuon::GetDS3()
{
  // project the road to station 3
  double Xproj = double(X[GP0] + Dx * (Z[ST3] - Z[GP0]) );
  double Yproj = double(Y[GP0] + Dy * (Z[ST3] - Z[GP0]) );

  // distance between projected position and station 3 point.
  double ds3 = sqrt( (Xproj-X[ST3])*(Xproj-X[ST3]) +
		     (Yproj-Y[ST3])*(Yproj-Y[ST3]) );
  return ds3;
}
//-----------------------------------------------------------------------------


inline
double IMuon::GetDDG0()
{
  // this is the angle between two vectors:  cos(t) = (A . B)/(|A||B|)
  // if the kalman projection is available the projected (Px, Py, Pz) is used.
  // vector at st3 is (Px, Py, Pz), vector at gap0 is (Dx, Dy, 1)*sign(Pz)

  int k=-9;
  // this is adabpted from MWG code
  // if kalman projection is not avaiavle them use momentum at st3.
  if (X[KAL] == 0 && Y[KAL] == 0 && Z[KAL] == 0) k = ST3;
  else k = KAL;

  double dz = (Pz[k] > 0) ? 1. : -1.; // select sign according to sign of Pz
  double dx = Dx*dz;   // this takes care of the sign fo Dx
  double dy = Dy*dz;

  double scaler = (Px[k]*dx + Py[k]*dy + Pz[k]*dz);
  scaler /= GetP(k);
  scaler /= sqrt(dx*dx + dy*dy + 1.); // dz*dz = 1

  if (scaler > 1.) scaler = 1.;

  return 57.295779506*acos(scaler); // 180/Pi = 57.3
}
//-----------------------------------------------------------------------------


inline
double IMuon::GetDS3cpt()
{
  return sqrt(pow((X[ST3] - X[GP0]* Z[ST3]/Z[GP0]), 2) +
	      pow((Y[ST3] - Y[GP0]* Z[ST3]/Z[GP0]), 2));
}
//-----------------------------------------------------------------------------


inline
short IMuon::GetNMuTrHits()
{
  short n=0;
  for (short j=15; j>=0; j--)
    if (short((MuTrHitPat >> j) & 1) == 1) n++;
  return n;
}
//-----------------------------------------------------------------------------


inline
short IMuon::GetNMuIdHits()
{
  short n=0;
  for (short j=9; j>=0; j--)
    if (short((MuIdHitPat >> j) & 1) == 1) n++;
  return n;
}
//-----------------------------------------------------------------------------


inline
void IMuon::Copy(IMuon *muon)
{
  for (short i=0; i<nPOSITIONS; i++)
    {
      X[i] = muon->X[i];
      Y[i] = muon->Y[i];
      Z[i] = muon->Z[i];
    }

  for (short i=0; i<nPOSITIONS-1; i++)
    {
      Px[i] = muon->Px[i];
      Py[i] = muon->Py[i];
      Pz[i] = muon->Pz[i];
    }

  Dx = muon->Dx;
  Dy = muon->Dy;

  Charge = muon->Charge;
  MuTrHitPat = muon->MuTrHitPat;
  MuIdHitPat = muon->MuIdHitPat;


  MuTrChi2 = muon->MuTrChi2;
  MuIdChi2 = muon->MuIdChi2;
  VtxRefitChi2 = muon->VtxRefitChi2;
}
//-----------------------------------------------------------------------------
#endif
