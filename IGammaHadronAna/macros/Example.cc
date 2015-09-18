
#include <iostream>
#include <fstream>
#include <assert.h>

#include "IPdst.hh"
#include "IConsts.hh"
#include "ITrack.hh"


#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"


using namespace std;

void Example()
{
  cout << "this is how to run IPdst analysis" << endl;

  assert(IConsts::Defined);
  IConsts::Print();

  ifstream fin(IConsts::RUNLIST);

  int runNum;

  TH1D *hPTA = new TH1D("hPTA","hPTA",100,IConsts::MINPTA,IConsts::MAXPTA);
  TH1D *hPC3sig = new TH1D("hPC3sig","hPC3sig",50,0,IConsts::MATCHSIG);
  TH1D *hEMCsig = new TH1D("hEMCsig","hEMCsig",50,0,IConsts::MATCHSIG);
  TH2D *hMatchSig = new TH2D("hMatchSig","hMatchSig",
			     50,-IConsts::MATCHSIG,IConsts::MATCHSIG,
			     50,-IConsts::MATCHSIG,IConsts::MATCHSIG);


  double pc3match, emcmatch;

  while (fin >> runNum)
    {
      IPdst *pdst = new IPdst(runNum, IConsts::INDIR);

      // totol number of events in a run.
      // IPdst joins all segments of a run.
      long N = pdst->GetEntries();
      cout << "run " << runNum << ",  entries " << N << endl;

      for (long i=0; i<N; i++)
	{
	  if (i%100000 == 0) cout << i << endl;
	  i++;
	  i++;

	  pdst->GetEntry(i);

	  int ntrks = pdst->_nTracks;  // number of track in event

	  for (int j=0; j<ntrks; j++)
	    {
	      ITrack *trk = pdst->_track[j]; // returns a pointer to ith track

	      // all the function of ITrack class can be accessed now.
	      if (fabs(trk->GetZedDC()) < 5.) continue;

	      if (trk->GetN0() > 0) continue;
 
	      pc3match = trk->GetPc3MatchingSig();
	      emcmatch = trk->GetEmcMatchingSig();

	      // PC3 and EMC matching
	      if (pc3match > IConsts::MATCHSIG && pc3match < 1e6) continue;
	      if (pc3match > 1e6 && emcmatch > IConsts::MATCHSIG) continue;

	      hPC3sig->Fill(pc3match);
	      hEMCsig->Fill(emcmatch);

	      hMatchSig->Fill(trk->GetPc3sdPhi(), trk->GetPc3sdZ());

	      hPTA->Fill(trk->GetPt());
	    }
	}

      cout << "---------------------------------------\n" << endl;
      delete pdst;
    }


  TCanvas *cc = new TCanvas("cc","cc",1000,800);
  cc->Divide(2,2);

  cc->cd(1);
  hPTA->Draw();

  cc->cd(2);
  hPC3sig->Draw();

  cc->cd(3);
  hEMCsig->Draw();

  cc->cd(4);
  hMatchSig->Draw("colz");
}
//-----------------------------------------------------------------------------
