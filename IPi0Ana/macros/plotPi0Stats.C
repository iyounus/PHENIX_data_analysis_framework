
{

  //ifstream fin("taxi103_runList_AN721short_AN689_prob01.txt");
  ifstream fin("../../ITextFiles/Run8pp_taxi174_runlist.txt");
  float run[1000];
  int nruns=0;

  while (fin >> run[nruns]) nruns++;

  cout << nruns << endl;

  fin.close();

  int runNum=0;
  float mean=0, meane=0, sigma=0, sigmae=0, sb=0;

  TFile *f = new TFile("Run8pp_ERT_pi0Stats_2Oct09_PbSc.root");
  TTree *t = (TTree *)f->Get("pi0");

  t->SetBranchAddress("runNum", &runNum);
  t->SetBranchAddress("mean",   &mean);
  t->SetBranchAddress("meane",  &meane);
  t->SetBranchAddress("sigma",  &sigma);
  t->SetBranchAddress("sigmae", &sigmae);
  t->SetBranchAddress("SB",     &sb);

  float Mean[1000];
  float MeanE[1000];
  float Sig[1000];
  float SigE[1000];
  float SB[1000];
  int index=0;


  for (int i=0; i<t->GetEntries(); i++)
    {
      t->GetEntry(i);

      //cout << runNum << endl;

      for (int j=0; j<nruns; j++)
	if ((int)run[j] == runNum)
 	  {
 	    Mean[j]  = mean;
 	    MeanE[j] = meane;
 	    Sig[j]   = sigma;
 	    SigE[j]  = sigmae;
	    SB[j]    = sb;
// 	    //if (sigma < 0.0082) cout << run[j] << "\t" << sigma << endl;
//  	    break;
	  }
    }

  TCanvas *cc =  new TCanvas("cc","cc",1200,600);
  cc->Divide(2,1);


  TGraphErrors *gMean = new TGraphErrors(nruns, run, Mean, 0, MeanE);
  TGraphErrors *gSig  = new TGraphErrors(nruns, run, Sig, 0, SigE);

  cc->cd(1);
  gMean->SetTitle("pi0 mass (PbSc);run number;mean");
  gMean->SetLineColor(41);
  gMean->SetMaximum(0.1395);
  gMean->SetMinimum(0.1365);
  gMean->GetXaxis()->CenterTitle();
  gMean->GetYaxis()->CenterTitle();
  gMean->Draw("AP");
  gMean->Draw("Psame");
  gMean->Fit("pol0");

  cc->cd(2);
  gSig->SetTitle("pi0 widths (PbSc);run number;sigma");
  gSig->SetMaximum(0.011);
  gSig->SetMinimum(0.009);
  gSig->SetLineColor(41);
  gSig->GetXaxis()->CenterTitle();
  gSig->GetYaxis()->CenterTitle();
  gSig->Draw("AP");
  gSig->Draw("Psame");
  gSig->Fit("pol0");

}

