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


void ChargeRatio_z()
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

  TFile *f = new TFile("IGammaHadron_pp.root");
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
  //double za[] = {0.,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1};
  const int nbins = 5;
  double za1[] = {0., 0.2, 0.4, 0.6, 0.8, 1.0};

  TString name = "";

  TH1D *hza[nbins];
  TH1D *hza_pos[nbins];
  TH1D *hza_neg[nbins];

  for (int i=0; i<nbins; i++)
    {
      name = "za";
      name += i;
      hza[i] = new TH1D(name.Data(), name.Data(), 50, za1[i], za1[i+1]);

      name = "za_pos";
      name += i;
      hza_pos[i] = new TH1D(name.Data(), name.Data(), 50, za1[i], za1[i+1]);

      name = "za_neg";
      name += i;
      hza_neg[i] = new TH1D(name.Data(), name.Data(), 50, za1[i], za1[i+1]);
    }

  TH1D *hmass = new TH1D("hmass","hmass",50,0,0.25);

  for (int i=0; i<gam->GetEntries(); i++)
    {
      if (i%100000 == 0) cout << i << endl;
      gam->GetEntry(i);

      if (phoPt < 4.5) continue;
      if (za < za1[0] || za > za1[nbins]) continue;
      //if (0.1*phoEn < eCone) continue;
      //if (prob > 0.2) continue;
      //if (diPhoMs > 0.097647 && diPhoMs < 0.17069) continue;


      int zabin = -9;
      for (int j=0; j<nbins; j++)
      	if (za > za1[j] && za < za1[j+1])
      	  {
      	    zabin = j;
      	    break;
      	  }
      if (zabin < 0) continue;


      hza[zabin]->Fill(za);

      if (hadCh > 0)
	hza_pos[zabin]->Fill(za);
      if (hadCh < 0)
	hza_neg[zabin]->Fill(za);

      hmass->Fill(diPhoMs);
    }


  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  c1->Divide(5,2);


  double mean_za[nbins];
  double ratio[nbins];
  double ratio_e[nbins];


  for (int i=0; i<nbins; i++)
    {
      mean_za[i] = hza[i]->GetMean();

      c1->cd(i+1);
      hza_pos[i]->Draw();
      double n1 = double(hza_pos[i]->GetEntries());

      c1->cd(i+nbins+1);
      hza_neg[i]->Draw();
      double n2 = double(hza_neg[i]->GetEntries());

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


  TGraphErrors *gza = new TGraphErrors(nbins, mean_za, ratio, 0, ratio_e);

  TCanvas *cc = new TCanvas("cc","cc",800,600);
  cc->cd();
  gza->SetTitle(";z_{F};R(h+/h-)");
  gza->SetMinimum(0.7);
  gza->SetMaximum(1.5);
  gza->GetXaxis()->SetLimits(0.,1.0);
  gza->SetMarkerColor(2);
  gza->SetMarkerStyle(22);
  gza->Draw("AP");

  TLine *l1 = new TLine(0., 1., 1., 1.);
  l1->SetLineStyle(2);
  l1->Draw();

}
