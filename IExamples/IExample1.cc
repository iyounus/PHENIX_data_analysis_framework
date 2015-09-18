/*
  June 20, 2010.
  This is a simple example of reading unm_pdst. In this example, I'm not
  creating any class. This is just one function.
*/

#include <iostream>
#include <assert.h>

#include "IPdst.hh"
#include "IConsts.hh"
#include "IHeader.hh"
#include "IHeaderList.hh"
#include "IPhoton.hh"
#include "IPhotonList.hh"

#include "TH1D.h"
#include "TLorentzVector.h"


using namespace std;

void IExample1(int runNum)
{
  /*
    IConsts must be definded before calling this function. This is the basic
    scheme for unm_pdst analysis. IConsts defines input/output directories,
    pt bins (if used), and several other consts which may be used in the
    analysis. The advantage of this scheme is that one can run several
    instances of the same code on a multicore machine without recompiling.
  */
  assert(IConsts::Defined); // If IConsts is not defined, this will exit.
  IConsts::Print(); // This prints out all the definition in IConsts

  // define histograms
  TH1D *hDiPhoton   = new TH1D("hDiPhoton","hDiPhoton", 100, 0.05, 0.25);

  TLorentzVector *vDiPhoton = new TLorentzVector(0.,0.,0.,0.);

  // now open the pdst
  IPdst *pdst = new IPdst(runNum); // this constructor takes one run

  const int Nentries = pdst->GetEntries();
  cout << "number of entries " << Nentries << endl;

  double P1, En1, Px1, Py1, Pz1;
  double P2, En2, Px2, Py2, Pz2;
  double mass;

  // now loop over all entries in the pdst
  for (int evt=0; evt<Nentries; evt++)
    {
      if (evt%100000 == 0) cout << evt << endl;
      pdst->GetEntry(evt); // just like you would do with a TTree

      // Now you can read all objects from the pdst.
      IHeader *hdr = pdst->GetHeader();
      if (hdr->GetZVertex() > 30.) continue;

      // pdst knows how many phtons/muons/tracks are present in an event.
      short nPhotons = pdst->GetNPhotons();
      if (nPhotons < 2) continue;

      // fill di photon mass histogram
      for (int i=0; i<nPhotons-1; i++)
	{
	  // now you can get all the photons
	  // see IPhoton class for the all available functions.
	  IPhoton *ph1 = pdst->GetPhoton(i);
	  if (!ph1->Fired4x4c()) continue;

	  // NOTE: all the values obtained from nDSTs are floats.
	  // all values calculated in IPhoton class are returned as double
	  En1 = double(ph1->GetEnergy());
	  Px1 = ph1->GetPx();
	  Py1 = ph1->GetPy();
	  Pz1 = ph1->GetPz();

	  for (int j=i; j<nPhotons; j++)
	    {
	      IPhoton *ph2 = pdst->GetPhoton(j);
	      if (!ph2->Fired4x4c()) continue;

	      En2 = double(ph2->GetEnergy());
	      Px2 = ph2->GetPx();
	      Py2 = ph2->GetPy();
	      Pz2 = ph2->GetPz();

	      double Easym = fabs( (En1 - En2)/(En1 + En2) );
	      if (Easym > 0.8) continue;
	      // Alternatively one can use IConsts::EAsym defined in IConsts

	      vDiPhoton->SetPxPyPzE(Px1+Px2, Py1+Py2, Pz1+Pz2, En1+En2);

	      mass =  vDiPhoton->M();
	      if (mass <= 0.05 || mass > 0.25) continue;

	      hDiPhoton->Fill(mass);
	    }
	}
    }

  // now you can delete pdst;
  delete pdst;

  hDiPhoton->Draw();
}
