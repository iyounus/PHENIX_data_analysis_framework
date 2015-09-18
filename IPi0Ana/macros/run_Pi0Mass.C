
#include "Pi0Mass.cc"

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"


void run_Pi0Mass()
{
  int runNum=0;
  float nPi0s=0, SS=0, BB=0, SB=0, chi2=0,
    mean=0.137, meane=0, sigma=0.09, sigmae=0, entries=0;

  //TFile *f = new TFile("junk.root", "recreate");
  TFile *f = new TFile("Run8pp_ERT_pi0Stats_25Feb10_PbGl.root", "recreate");
  TTree *t = new TTree("pi0", "pi0");


  t->Branch("runNum", &runNum, "runNum/I");
  t->Branch("nPi0s",  &nPi0s,  "nPi0s/F");
  t->Branch("SS",     &SS,     "SS/F");
  t->Branch("BB",     &BB,     "BB/F");
  t->Branch("SB",     &SB,     "SB/F");
  t->Branch("chi2",   &chi2,   "chi2/F");
  t->Branch("mean",   &mean,   "mean/F");
  t->Branch("sigma",  &sigma,  "sigma/F");
  t->Branch("meane",  &meane,  "meane/F");
  t->Branch("sigmae", &sigmae, "sigmae/F");
  t->Branch("entires",&entries,"entries/F");

  TF1 *fn = 0;

  ifstream fin("../../ITextFiles/taxi191_list.txt");
  while (fin >> runNum)
    {
      mean = 0.135;
      sigma = 0.08;

      TString file = "/data1/IPdstPi0Histos_Run8pp200ERT/IPi0Mass_";
      file += runNum;
      file += ".root";

      cout << file.Data() << endl;

      Pi0Mass *ms1 = new Pi0Mass("1", "PbGl", file.Data());
      entries = float(ms1->Get_hMass()->GetEntries());
      cout << entries << endl;
      do
  	{
	  if (fn) delete fn;
	  fn = ms1->Fit_hMass((char*)("QRENO"), mean, sigma);

	  mean   = float(fn->GetParameter(1));
	  sigma  = float(fabs(fn->GetParameter(2)));
	  meane  = float(fn->GetParError(1));
	  sigmae = float(fn->GetParError(2));
	  chi2   = float(fn->GetChisquare()/double(fn->GetNDF()));

	  if (chi2 > 3.0) cout << runNum << endl;
   	}
      while (chi2 > 4.0);

      nPi0s = float(ms1->Get_nPi0s(9));

      SS = float(ms1->Get_SS(9));
      BB = float(ms1->Get_BB(9));
      SB = float(ms1->Get_SB(9));

      cout << mean << "\t" << sigma << "\t" 
	   <<  SS << "\t" << BB << "\t" << SB << endl;

      t->Fill();
      delete ms1;
    }


  f->Write();
  f->Close();
  delete f;
}
