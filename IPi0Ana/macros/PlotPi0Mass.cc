
void PlotPi0Mass(int runNum)
{
  gROOT->Reset();

  gSystem->Load("libCore.so");
  gSystem->Load("libCint.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");

  gSystem->AddIncludePath(" -I/home/iyounus/workdir/analysis/unm_pdst_ana/IContainers");

  gSystem->Load("libIContainers.so");
  gROOT->ProcessLine(".L Pi0Mass.cc+");

  TString file = "histograms_taxi103_ERT/pi0mass_";
  file += runNum;
  file += ".root";


  Pi0Mass *ms1 = new Pi0Mass("1", "PbSc", file.Data());

  float mean=0.135, sigma=0.08;
  TF1 *fn = ms1->Fit_hMass("R", mean, sigma);

}
