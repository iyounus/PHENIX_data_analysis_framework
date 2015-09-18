


void PlotPi0Stats()
{
  TFile *f = new TFile("Pi0Stats_02Apr06.root");
  TTree *t = (TTree*)f->Get("pi0");

  float mean, rms, nPi0, nBkg;
  float meanF, sigF, nPi0F, meanE, sigE, nPi0E, chi2;
  float par[5], err[5];
  int run, PtBin, PbSc;

  float minPt[3] = {1,3,6};
  float maxPt[3] = {3,6,15};
  TString det[2] = {"PbGl","PbSc"};

  t->SetBranchAddress("nPi0", &nPi0);
  t->SetBranchAddress("nBkg", &nBkg);
  t->SetBranchAddress("Mean", &mean);
  t->SetBranchAddress("RMS",  &rms);
  t->SetBranchAddress("MeanF",&meanF);
  t->SetBranchAddress("SigF", &sigF);
  t->SetBranchAddress("nPi0F",&nPi0F);
  t->SetBranchAddress("MeanE",&meanE);
  t->SetBranchAddress("SigE", &sigE);
  t->SetBranchAddress("nPi0E",&nPi0E);
  t->SetBranchAddress("Chi2", &chi2);
  t->SetBranchAddress("Run",  &run);
  t->SetBranchAddress("PtBin",&PtBin);
  t->SetBranchAddress("PbSc", &PbSc);

  float varX[10000], varY[10000], varXe[10000], varYe[10000];

  ifstream fin2("runList_short.txt");
  int run2[1000];
  int count2 = 0;
  while (fin2 >> run2[count2])
    count2++;

  fin2.close();

  cout << count2 << endl;

  ofstream fout("allRuns_Pi0Stats.txt");

  int N=0;
  for (int i=0; i<t->GetEntries(); i++)
    {
      t->GetEntry(i);
      bool include  = false;

      for (int j=0; j<count2; j++)
	if (run == run2[j]) include = true;

      float del=0;
      if (PbSc == 0 || PtBin != 1) continue;

      if ((nPi0+nBkg) != 0)
      	del = rms/TMath::Sqrt(nPi0+nBkg);

      if (include && mean > 0.12 && mean < 0.15 && del > 0 && del < 0.0025)
	{
	  //cout << run << "\t" << meanF << endl;
	  varX[N] = float(run);
	  varXe[N] = 0;

	  varY[N] = mean*1000;
	  varYe[N] = del*1000;

	  fout << run << "\t" << mean << "\t" << rms << endl;
	  N++;
	}
      else
	cout << run << "\t" << mean << "\t" << rms << endl;
    }
  //fout.close();
  cout << N << endl;
  TGraphErrors *gr = new TGraphErrors(N, varX, varY, varXe, varYe);
  gr->SetTitle("Pi0 mass (MeV), PbSc, 3 < pT < 6; run number");
  //gr->SetMaximum(139);
  //gr->SetMinimum(134);
  gr->Draw("AP");
  gr->Draw("P");


}
