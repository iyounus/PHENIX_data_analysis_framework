#include <iostream>

#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TLine.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"


void ChargeRatio_xE()
{
  //gStyle->SetStatW(0.5);

  double phoPhi;
  double hadPhi;
  double phoPt;
  double hadPt;
  double phoEn;
  double EP;
  double za;
  double xE;
  double dPhi;
  double eCone;
  double hadCh;
  double diHadMs;
  double diPhoMs;
  double prob;

  TFile *f = new TFile("IGammaHadron_dAu.root");
  TTree *gam = (TTree *)f->Get("gam");

  gam->SetBranchAddress("PhotonPt",  &phoPt);
  gam->SetBranchAddress("HadronPt",  &hadPt);
  gam->SetBranchAddress("PhotonPhi", &phoPhi);
  gam->SetBranchAddress("HadronPhi", &hadPhi);
  gam->SetBranchAddress("PhotonEn",  &phoEn);
  gam->SetBranchAddress("EP",        &EP);
  gam->SetBranchAddress("za",        &za);
  gam->SetBranchAddress("xE",        &xE);
  gam->SetBranchAddress("delPhi",    &dPhi);
  gam->SetBranchAddress("ECone",     &eCone);
  gam->SetBranchAddress("HadronCh",  &hadCh);
  gam->SetBranchAddress("DiHadMs",   &diHadMs);
  gam->SetBranchAddress("DiPhoMs",   &diPhoMs);
  gam->SetBranchAddress("Prob",      &prob);


  //const int nbins = 10;
  //double xe[] = {0.,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1};
  const int nbins = 5;
  double xe[] = {0., 0.2, 0.4, 0.6, 0.8, 1.0};

  TString name = "";

  TH1D *hxE[nbins];
  TH1D *hxE_pos[nbins];
  TH1D *hxE_neg[nbins];

  for (int i=0; i<nbins; i++)
    {
      name = "xE";
      name += i;
      hxE[i] = new TH1D(name.Data(), name.Data(), 50, xe[i], xe[i+1]);

      name = "xE_pos";
      name += i;
      hxE_pos[i] = new TH1D(name.Data(), name.Data(), 50, xe[i], xe[i+1]);

      name = "xE_neg";
      name += i;
      hxE_neg[i] = new TH1D(name.Data(), name.Data(), 50, xe[i], xe[i+1]);
    }

  TH1D *hmass = new TH1D("hmass","hmass",50,0,0.25);

  for (int i=0; i<gam->GetEntries(); i++)
    {
      if (i%100000 == 0) cout << i << endl;
      gam->GetEntry(i);

      if (phoPt < 4.) continue;
      if (xE < xe[0] || xE > xe[nbins]) continue;
      //if (0.1*phoEn < eCone) continue;
      if (prob > 0.2) continue;
      //if (diPhoMs > 0.097647 && diPhoMs < 0.17069) continue;


      int xebin = -9;
      for (int j=0; j<nbins; j++)
      	if (xE > xe[j] && xE < xe[j+1])
      	  {
      	    xebin = j;
      	    break;
      	  }
      if (xebin < 0) continue;


      hxE[xebin]->Fill(xE);

      if (hadCh > 0)
	hxE_pos[xebin]->Fill(xE);
      if (hadCh < 0)
	hxE_neg[xebin]->Fill(xE);

      hmass->Fill(diPhoMs);
    }


  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  c1->Divide(5,2);


  double mean_xE[nbins];
  double ratio[nbins];
  double ratio_e[nbins];


  for (int i=0; i<nbins; i++)
    {
      mean_xE[i] = hxE[i]->GetMean();

      c1->cd(i+1);
      hxE_pos[i]->Draw();
      double n1 = double(hxE_pos[i]->GetEntries());

      c1->cd(i+nbins+1);
      hxE_neg[i]->Draw();
      double n2 = double(hxE_neg[i]->GetEntries());

      ratio[i] = n1/n2;
      ratio_e[i] = ratio[i]*sqrt(1/n1 + 1/n2);
    }

  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->cd();
  hmass->SetTitle("di-photon mass;M_{#gamma#gamma} GeV");
  hmass->Draw();

  TF1 *fM = new TF1("fM","gaus(0) + pol1(3)",0.06,0.21);
  fM->SetParameter(0, 100.);
  fM->SetParameter(1, 0.135);
  fM->SetParameter(2, 0.02);
  fM->SetParameter(3, 10.);
  fM->SetParameter(4, 0.001);

  hmass->Fit(fM, "R");
  double meanM = fM->GetParameter(1);
  double sigma = fM->GetParameter(2);
  double npi0  = fM->Integral(meanM - 3*sigma, meanM + 3*sigma) /
    hmass->GetBinWidth(10);
  cout << npi0 << endl;
  cout << meanM - 3*sigma << "  " << meanM + 3*sigma << endl;


  TGraphErrors *gxE = new TGraphErrors(nbins, mean_xE, ratio, 0, ratio_e);

  TCanvas *cc = new TCanvas("cc","cc",800,600);
  cc->cd();
  gxE->SetTitle(";z_{F};R(h+/h-)");
  gxE->SetMinimum(0.7);
  gxE->SetMaximum(1.5);
  gxE->GetXaxis()->SetLimits(0.,1.0);
  gxE->SetMarkerColor(2);
  gxE->SetMarkerStyle(22);
  gxE->Draw("AP");

  TLine *l1 = new TLine(0., 1., 1., 1.);
  l1->SetLineStyle(2);
  l1->Draw();

}
