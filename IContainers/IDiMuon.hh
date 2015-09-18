//-----------------------------------------------------------
// IDiMuon.hh
//	 Created  Feb 23 2011
//       Imran Younus
//-----------------------------------------------------------

#ifndef IDIMUON_H
#define IDIMUON_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <math.h>

#define nPOSITIONS 6 // vtx, st1, st2, st3, kalman projection to gap0, gap0

#define MuMass   0.105658367     // from PDG
#define MuMass2  0.0111636905171 // square of muon mass


class IDiMuon : public TObject
{
private:

  float X[2][nPOSITIONS]; // vtx, st1, st2, st3, kalman proj, gap0
  float Y[2][nPOSITIONS];
  float Z[2][nPOSITIONS];

  float Px[2][nPOSITIONS-1]; // no momentum at gap0
  float Py[2][nPOSITIONS-1];
  float Pz[2][nPOSITIONS-1];

  float Dx[2];  // dxdz at MuId gap0
  float Dy[2];  // dydz at MuId gap0

  short Charge[2];

  unsigned short MuTrHitPat[2]; // hit pattern in MuTr
  unsigned short MuIdHitPat[2]; // hit pattern in MuId

  float MuTrChi2[2];
  float MuIdChi2[2];

  float VtxRefitChi2;
  float VtxRefitZ;

public:
  IDiMuon(){;}
  ~IDiMuon(){;}

  enum POS  {VTX, ST1, ST2, ST3, KAL, GP0};
  enum MUON {Muon1, Muon2};
  enum MU   {Mu1, Mu2};  // I don't like to write Muon all the time!

  float GetX(short imu, short pos)   const {return X[imu][pos];}
  float GetY(short imu, short pos)   const {return Y[imu][pos];}
  float GetZ(short imu, short pos)   const {return Z[imu][pos];}

  float GetPx(short imu, short pos)  const {return Px[imu][pos];}
  float GetPy(short imu, short pos)  const {return Py[imu][pos];}
  float GetPz(short imu, short pos)  const {return Pz[imu][pos];}

  float GetDx(short imu)             const {return Dx[imu];};
  float GetDy(short imu)             const {return Dy[imu];};
  float GetMuIdSlopeX(short imu)     const {return Dx[imu];};
  float GetMuIdSlopeY(short imu)     const {return Dy[imu];};

  short GetCharge(short imu)         const {return Charge[imu];}

  unsigned short GetMuTrHitPat(short imu) const {return MuTrHitPat[imu];}
  unsigned short GetMuIdHitPat(short imu) const {return MuIdHitPat[imu];}

  float GetMuTrChi2(short imu)       const {return MuTrChi2[imu];}
  float GetMuIdChi2(short imu)       const {return MuIdChi2[imu];}

  float GetVtxRefitChi2()            const {return VtxRefitChi2;}
  float GetVtxRefitZ()               const {return VtxRefitZ;}


  // dimuon function
  double GetPx()    const {return double(Px[Muon1][VTX] + Px[Muon2][VTX]);}
  double GetPy()    const {return double(Py[Muon1][VTX] + Py[Muon2][VTX]);}
  double GetPz()    const {return double(Pz[Muon1][VTX] + Pz[Muon2][VTX]);}
  double GetP2()    // momentum square
  {return GetPx()*GetPx() + GetPy()*GetPy() + GetPz()*GetPz();}

  double GetMass();
  double GetE()           {return GetE(Muon1) + GetE(Muon2);}
  double GetP()           {return sqrt(GetP2());}
  double GetPt()          {return sqrt(GetPx()*GetPx() + GetPy()*GetPy());}
  double GetPhi()         {return atan2(GetPy(), GetPx());}
  double GetTheta()       {return atan2(GetPt(), GetPz());}
  double GetRapidity()
  {return 0.5*log( (GetE() + GetPz())/(GetE() - GetPz()) );}

  // these funstion take center-of-mass energy as input (e.g., ss = 200)
  double GetX1(double ss) {return (GetMass()/ss) * exp( GetRapidity()); }
  double GetX2(double ss) {return (GetMass()/ss) * exp(-GetRapidity()); }
  double GetXF(double ss) {return GetX1(ss) - GetX2(ss);}



  // single muon functions
  double GetE        (short imu);
  double GetP        (short imu, short pos);
  double GetPt       (short imu, short pos);

  double GetPhi      (short imu, short pos)
  {return atan2(Py[imu][pos], Px[imu][pos]);}

  double GetTheta    (short imu, short pos)
  {return atan2(GetPt(imu, pos), Pz[imu][pos]);}

  // Eta and Rapidity are the same, just giving a name choice to user
  double GetEta      (short imu);
  double GetRapidity (short imu);
  double GetELoss    (short imu)
  {return GetP(imu, VTX) - GetP(imu, ST1);}

  short  GetNMuTrHits(short imu);
  short  GetNMuIdHits(short imu);
  double GetDG0      (short imu);
  double GetDDG0     (short imu);
  double GetDS3      (short imu);
  double GetDS3cpt   (short imu);



  // single muon function
  double GetE_Muon1()              {return GetE(Muon1);}
  double GetE_Muon2()              {return GetE(Muon2);}

  double GetP_Muon1  (short pos)   {return GetP(Muon1, pos);}
  double GetP_Muon2  (short pos)   {return GetP(Muon2, pos);}

  double GetPt_Muon1 (short pos)   {return GetPt(Muon1, pos);}
  double GetPt_Muon2 (short pos)   {return GetPt(Muon2, pos);}

  double GetPhi_Muon1(short pos)   {return GetPhi(Muon1, pos);}
  double GetPhi_Muon2(short pos)   {return GetPhi(Muon2, pos);}

  double GetTheta_Muon1(short pos) {return GetTheta(Muon1, pos);}
  double GetTheta_Muon2(short pos) {return GetTheta(Muon2, pos);}

  double GetELoss_Muon1()          {return GetELoss(Muon1);}
  double GetELoss_Muon2()          {return GetELoss(Muon2);}

  double GetEta_Muon1()            {return GetEta(Muon1);}
  double GetEta_Muon2()            {return GetEta(Muon2);}
  double GetRapidity_Muon1()       {return GetRapidity(Muon1);}
  double GetRapidity_Muon2()       {return GetRapidity(Muon2);}

  short  GetNMuTrHits_Muon1()      {return GetNMuTrHits(Muon1);}
  short  GetNMuTrHits_Muon2()      {return GetNMuTrHits(Muon2);}
  short  GetNMuIdHits_Muon1()      {return GetNMuIdHits(Muon1);}
  short  GetNMuIdHits_Muon2()      {return GetNMuIdHits(Muon2);}

  double GetDG0_Muon1()            {return GetDG0(Muon1);}
  double GetDG0_Muon2()            {return GetDG0(Muon2);}
  double GetDDG0_Muon1()           {return GetDDG0(Muon1);}
  double GetDDG0_Muon2()           {return GetDDG0(Muon2);}
  double GetDS3_Muon1()            {return GetDS3(Muon1);}
  double GetDS3_Muon2()            {return GetDS3(Muon2);}
  double GetDS3cpt_Muon1()         {return GetDS3cpt(Muon1);}
  double GetDS3cpt_Muon2()         {return GetDS3cpt(Muon2);}



  void SetX(short i, short pos, float x)   {X[i][pos] = x;}
  void SetY(short i, short pos, float y)   {Y[i][pos] = y;}
  void SetZ(short i, short pos, float z)   {Z[i][pos] = z;}

  void SetPx(short i, short pos, float px) {Px[i][pos] = px;}
  void SetPy(short i, short pos, float py) {Py[i][pos] = py;}
  void SetPz(short i, short pos, float pz) {Pz[i][pos] = pz;}

  void SetDx(short i, float dx) {Dx[i] = dx;}
  void SetDy(short i, float dy) {Dy[i] = dy;}

  void SetCharge(short i, short ch)    {Charge[i] = ch;}

  void SetMuTrHitPat(short i, unsigned short n) {MuTrHitPat[i] = n;}
  void SetMuIdHitPat(short i, unsigned short n) {MuIdHitPat[i] = n;}


  void SetMuTrChi2(short i, float chi) {MuTrChi2[i] = chi;}
  void SetMuIdChi2(short i, float chi) {MuIdChi2[i] = chi;}

  void SetVtxRefitChi2(float chi2)     {VtxRefitChi2 = chi2;}
  void SetVtxRefitZ   (float z)        {VtxRefitZ = z;}

  void Copy(IDiMuon *muon);

protected:
  ClassDef(IDiMuon, 2)
};
//=============================================================================


// inline functions

inline 
double IDiMuon::GetMass()
{
  double E1 = GetE(Muon1);
  double E2 = GetE(Muon2);
  double M  = sqrt( (E1+E2)*(E1+E2) - GetP2());
  return M;
}
//-----------------------------------------------------------------------------


inline
double IDiMuon::GetE(short imu)
{
  return sqrt( Px[imu][VTX]*Px[imu][VTX] +
	       Py[imu][VTX]*Py[imu][VTX] +
	       Pz[imu][VTX]*Pz[imu][VTX] +
	       MuMass2 );
}
//-----------------------------------------------------------------------------


inline
double IDiMuon::GetP(short imu, short pos)
{
  return sqrt( Px[imu][pos]*Px[imu][pos] +
	       Py[imu][pos]*Py[imu][pos] +
	       Pz[imu][pos]*Pz[imu][pos] );
}
//-----------------------------------------------------------------------------


inline
double IDiMuon::GetPt(short imu, short pos)
{
  return sqrt( Px[imu][pos]*Px[imu][pos] + Py[imu][pos]*Py[imu][pos] );
}
//-----------------------------------------------------------------------------


inline
double IDiMuon::GetRapidity(short imu)
{
  double eta = 0.5*log( (GetE(imu) + Pz[imu][VTX]) /
			(GetE(imu) - Pz[imu][VTX]) );
  return eta;
}
//-----------------------------------------------------------------------------


inline
double IDiMuon::GetEta(short imu)
{
  double eta = 0.5*log( (GetE(imu) + Pz[imu][VTX]) /
			(GetE(imu) - Pz[imu][VTX]) );
  return eta;
}
//-----------------------------------------------------------------------------


inline
double IDiMuon::GetDG0(short imu)
{
  double Xproj = 9.e12;
  double Yproj = 9.e12;

  if (X[imu][KAL] == 0 && Y[imu][KAL] == 0)
    {
      Xproj = double(X[imu][ST3] + Px[imu][ST3]*
		     (Z[imu][GP0] - Z[imu][ST3]) /Pz[imu][ST3]);

      Yproj = double(Y[imu][ST3] + Py[imu][ST3]*
		     (Z[imu][GP0] - Z[imu][ST3]) /Pz[imu][ST3]);
    }
  else
    {
      Xproj = double(X[imu][KAL] + Px[imu][ST3]*
		     (Z[imu][GP0] - Z[imu][KAL]) /Pz[imu][ST3]);

      Yproj = double(Y[imu][KAL] + Py[imu][ST3]*
		     (Z[imu][GP0] - Z[imu][KAL]) /Pz[imu][ST3]);
    }


  double dg0 = sqrt( (Xproj-X[imu][GP0])*(Xproj-X[imu][GP0]) +
		     (Yproj-Y[imu][GP0])*(Yproj-Y[imu][GP0]) );
  return dg0;
}
//-----------------------------------------------------------------------------


inline
double IDiMuon::GetDS3(short imu)
{
  // project the road to station 3
  double Xproj = double(X[imu][GP0] + Dx[imu] * (Z[imu][ST3] - Z[imu][GP0]));
  double Yproj = double(Y[imu][GP0] + Dy[imu] * (Z[imu][ST3] - Z[imu][GP0]));

  // distance between projected position and station 3 point.
  double ds3 = sqrt( (Xproj-X[imu][ST3])*(Xproj-X[imu][ST3]) +
		     (Yproj-Y[imu][ST3])*(Yproj-Y[imu][ST3]) );
  return ds3;
}
//-----------------------------------------------------------------------------


inline
double IDiMuon::GetDDG0(short imu)
{
  // this is the angle between two vectors:  cos(t) = (A . B)/(|A||B|)
  // if the kalman projection is available the projected (Px, Py, Pz) is used.
  // vector at st3 is (Px, Py, Pz), vector at gap0 is (Dx, Dy, 1)*sign(Pz)

  int k=-9;
  // this is adabpted from MWG code
  // if kalman projection is not avaiavle them use momentum at st3.
  if (X[imu][KAL] == 0 && Y[imu][KAL] == 0 && Z[imu][KAL] == 0) k = ST3;
  else k = KAL;

 // select sign according to sign of Pz
  double dz = (Pz[imu][k] > 0) ? 1. : -1.;
  double dx = Dx[imu]*dz;   // this takes care of the sign fo Dx
  double dy = Dy[imu]*dz;

  double scaler = (Px[imu][k]*dx + Py[imu][k]*dy + Pz[imu][k]*dz);
  scaler /= GetP(imu, k);
  scaler /= sqrt(dx*dx + dy*dy + 1.); // dz*dz = 1

  if (scaler > 1.) scaler = 1.;

  return 57.295779506*acos(scaler); // 180/Pi = 57.3
}
//-----------------------------------------------------------------------------


inline
double IDiMuon::GetDS3cpt(short imu)
{
  return sqrt(pow((X[imu][ST3] - X[imu][GP0]* Z[imu][ST3]/Z[imu][GP0]), 2) +
	      pow((Y[imu][ST3] - Y[imu][GP0]* Z[imu][ST3]/Z[imu][GP0]), 2));
}
//-----------------------------------------------------------------------------


inline
short IDiMuon::GetNMuTrHits(short imu)
{
  short n=0;
  for (short j=15; j>=0; j--)
    if (short((MuTrHitPat[imu] >> j) & 1) == 1) n++;
  return n;
}
//-----------------------------------------------------------------------------


inline
short IDiMuon::GetNMuIdHits(short imu)
{
  short n=0;
  for (short j=9; j>=0; j--)
    if (short((MuIdHitPat[imu] >> j) & 1) == 1) n++;
  return n;
}
//-----------------------------------------------------------------------------


inline
void IDiMuon::Copy(IDiMuon *dimuon)
{
  for (int muon=0; muon<2; muon++)
    {
      for (short pos=0; pos<nPOSITIONS;pos++)
	{
	  X[muon][pos] = dimuon->X[muon][pos];
	  Y[muon][pos] = dimuon->Y[muon][pos];
	  Z[muon][pos] = dimuon->Z[muon][pos];
	}

      for (short pos=0; pos<nPOSITIONS-1; pos++)
	{
	  Px[muon][pos] = dimuon->Px[muon][pos];
	  Py[muon][pos] = dimuon->Py[muon][pos];
	  Pz[muon][pos] = dimuon->Pz[muon][pos];
	}

      Dx[muon] = dimuon->Dx[muon];
      Dy[muon] = dimuon->Dy[muon];

      Charge[muon] = dimuon->Charge[muon];
      MuTrHitPat[muon] = dimuon->MuTrHitPat[muon];
      MuIdHitPat[muon] = dimuon->MuIdHitPat[muon];

      MuTrChi2[muon] = dimuon->MuTrChi2[muon];
      MuIdChi2[muon] = dimuon->MuIdChi2[muon];
    }

  VtxRefitChi2 = dimuon->VtxRefitChi2;
  VtxRefitZ    = dimuon->VtxRefitZ;
}
//-----------------------------------------------------------------------------
#endif
